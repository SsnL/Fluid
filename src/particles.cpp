#include <sstream>

#include "particles.h"
#include "misc/sphere_drawing.h"
#include "static_scene/marching_triangle.h"

// params
// particle visualization radius
#define VIS_RADIUS 0.05
// the threshold of isovalue for isosurface
#define ISO_LEVEL 700
// the resolution of isosurface
#define GRADIENT_EPS 0.001
// the resolution of isosurface
#define FSTEPSIZE 0.18
// velocity bounce factor
#define VELOCITY_BOUNCE_FACTOR 0.6
// default time step
#define DEFAULT_DELTA_T 0.016
// used in SPH density estimation
#define H 0.3
// h^2, neighbor radius squared
#define H2 0.09
// max num neighbors
// #define MAX_NUM_NEIGHBORS 45
// should be 30~40 for SPH density est to be stable
#define NUM_NEIGHBOR_ALERT_THRESHOLD 18
// number of Newton steps
#define NEWTON_NUM_STEPS 12
// used in calculating lambda
#define EPSILON 2
// artificial pressure coeff
#define K 0.0001
// artificial pressure exp
#define N 4
// vorticity
#define VORTICITY_EPSILON 0.001
// XSPH viscocity
#define C 0.001
// num cells along axis
// #define CELL_NUM_X 8
// #define CELL_NUM_Y 6
// #define CELL_NUM_Z 8

namespace CGL {
  inline void clamp(Vector3D &p, Vector3D delta_p, BVHAccel *bvh) {
    double l = delta_p.norm();
    if (l <= EPS_D) {
      return;
    }
    Vector3D d = delta_p / l;
    // Viewing from plane x = -1.0, clip it as well
    bool intersect = false;
    if (d.z != 0.0) {
      double pt = (1.0 - p.z) / d.z;
      if (pt > 0.0 && pt < l) {
        l = pt;
        intersect = true;
      }
    }
    // light is at y = 1.49.
    // to avoid particles being stuck between light and ceiling, clip it
    if (d.y != 0.0) {
      double pt = (1.49 - p.y) / d.y;
      if (pt > 0.0 && pt < l) {
        l = pt;
        intersect = true;
      }
    }
    Ray r(p, d, 0.0, l);
    if (bvh->intersect(r) || intersect) {
      p += (r.max_t  - EPS_D) * d;
    } else {
      p += delta_p;
    }
    p.x = max(-1.0 + EPS_D, min(1.0 - EPS_D, p.x));
    p.y = max(0.0 + EPS_D, min(1.49 - EPS_D, p.y));
    p.z = max(-1.0 + EPS_D, min(1.0 - EPS_D, p.z));
  }

  // Set non-zero velocity response if collision
  inline void clamp_response(Vector3D &p, Vector3D &v, double delta_t, BVHAccel *bvh) {
    Vector3D delta_p = v * delta_t;
    double total_l = delta_p.norm(), l = total_l;
    if (l <= EPS_D) {
      return;
    }
    Vector3D d = delta_p / l;
    // Viewing from plane x = -1.0, clip it as well
    bool intersect = false;
    if (d.z != 0.0) {
      double pt = (1.0 - p.z) / d.z;
      if (pt > 0.0 && pt < l) {
        l = pt;
        intersect = true;
      }
    }
    // light is at y = 1.49.
    // to avoid particles being stuck between light and ceiling, clip it
    if (d.y != 0.0) {
      double pt = (1.49 - p.y) / d.y;
      if (pt > 0.0 && pt < l) {
        l = pt;
        intersect = true;
      }
    }
    Ray r(p, d, 0.0, l);
    Intersection i;
    if (bvh->intersect(r, &i) || intersect) {
      p += (r.max_t - EPS_D) * d;
      // hack. redirect ONCE
      // if not perpendicular
      if (1 - dot(d, i.n) > EPS_D && !intersect) { //TODO: INTERSECT = true case
        i = Intersection();
        d = (delta_p - dot(delta_p, i.n) * i.n).unit();
        r = Ray(p, d, 0.0, (total_l - r.max_t));
        bvh->intersect(r, &i);
        p += (r.max_t - EPS_D) * d;
      }

    } else {
      p += delta_p;
    }
    p.x = max(-1.0 + EPS_D, min(1.0 - EPS_D, p.x));
    p.y = max(0.0 + EPS_D, min(1.49 - EPS_D, p.y));
    p.z = max(-1.0 + EPS_D, min(1.0 - EPS_D, p.z));
  }

  inline double poly6_kernel(Vector3D r) {
    double r_l_2 = r.norm2();
    if (r_l_2 >= H2) {
      return 0.0;
    }
    double temp = H2 - r_l_2;
    return 1.56668147106 * intpow<3>(temp) / intpow<9>(H);
  }

  inline Vector3D grad_spiky_kernel(Vector3D r) {
    double r_l = r.norm();
    if (r_l >= H || r_l < EPS_D) {
      return Vector3D();
    }
    return -3 * 4.774648292756860 * intpow<2>(H - r_l) * r / (intpow<6>(H) * r_l);
  }

  double Particle::tensile_instability_scale = 1.0 / poly6_kernel(Vector3D(0, 0, 0.1 * H));

  std::ostream& operator<<( std::ostream& os, Particle& v ) {
    os << "P(p" << v.getNewPosition() << ",v" << v.velocity << ')';
    return os;
  }

  void Particle::estimateDensityWithNeighbors(std::vector<Particle *> neighbors) {
    density = 0.0;
    for (Particle *n: neighbors) {
      density += poly6_kernel(n->position - position);
    }
  }

  void Particle::initializeWithNewNeighbors() {
    s_corr.resize(neighbors.size());
    grad_w_neighbors.resize(neighbors.size());
    if (neighbors.size() < NUM_NEIGHBOR_ALERT_THRESHOLD) {
      cerr << *this
        << " only has " << neighbors.size()
        << " neighbors." << endl;
    }
  }

  void Particle::applyForceVelocity(BVHAccel *bvh, double delta_t) {
    // for (Force *f : particles->fs) {
    //   velocity += delta_t * f->getAccerlation(particles->simulate_time, *this);
    // }
    // assume only gravity for now
    velocity.y -= 10 * delta_t;
    new_position = position;
    clamp_response(new_position, velocity, delta_t, bvh);
  }

  void Particle::newtonStepCalculateLambda() {
    density = 0.0;
    double denom = 0.0;
    Vector3D grad_i_c_i;
    for (size_t i = 0; i < neighbors.size(); i++) {
      Vector3D r = new_position - neighbors[i]->new_position;
      double w_i_n = poly6_kernel(r);
      Vector3D grad_w_i_n = grad_spiky_kernel(r);
      // SPH density esitmator
      density += w_i_n;
      s_corr[i] = -K * intpow<N>(w_i_n * tensile_instability_scale);
      grad_w_neighbors[i] = grad_w_i_n;
      grad_w_i_n /= rest_density;
      grad_i_c_i += grad_w_i_n;
      denom += grad_w_i_n.norm2();
    }
    double c_i = density / rest_density - 1.0;
    denom += grad_i_c_i.norm2();
    lambda = -c_i / (denom + EPSILON);
  }

  void Particle::newtonStepUpdatePosition(BVHAccel *bvh) {
    Vector3D delta_p;
    for (size_t i = 0; i < neighbors.size(); i++) {
      delta_p += (lambda + neighbors[i]->lambda + s_corr[i]) * grad_w_neighbors[i];
    }
    delta_p /= rest_density;
    clamp(new_position, delta_p, bvh);
  }

  void Particle::updateVelocity(double delta_t) {
    velocity = (new_position - position) / delta_t;
  }

  void Particle::calculateVorticityApplyXSPHViscosity() {
    vorticity = Vector3D();
    Vector3D viscocity;
    density = 0.0;
    for (size_t i = 0; i < neighbors.size(); i++) {
      Vector3D v_i_n = neighbors[i]->velocity - velocity;
      Vector3D r = new_position - neighbors[i]->new_position;
      Vector3D grad_w_i_n = grad_spiky_kernel(r);
      vorticity += cross(v_i_n, grad_w_i_n);
      double w_i_n = poly6_kernel(r);
      viscocity += v_i_n * w_i_n;
      grad_w_neighbors[i] = grad_w_i_n;
      density += w_i_n;
    }
    velocity += C * viscocity;
  }

  void Particle::applyVorticity(double delta_t) {
    Vector3D grad_vorticity;
    for (size_t i = 0; i < neighbors.size(); i++) {
      grad_vorticity += neighbors[i]->vorticity.norm() * grad_w_neighbors[i];
    }
    if (grad_vorticity.norm() > EPS_D) {
      velocity += delta_t * VORTICITY_EPSILON * cross(grad_vorticity.unit(), vorticity);
    }
  }

  void Particle::updatePosition() {
    position = new_position;
  }

  void Particles::timeStep(double delta_t) {
    cerr << "Time: " << simulate_time;
    simulate_time += delta_t;
    cerr << " => " << simulate_time << endl;
    for (Particle *p : ps) {
      p->neighbors.clear();
      p->applyForceVelocity(bvh, delta_t);
    }
    for (size_t i = 0; i < ps.size(); i++) {
      for (size_t j = i + 1; j < ps.size(); j++) {
        if ((ps[i]->getNewPosition() - ps[j]->getNewPosition()).norm2() <= H2) {
          ps[i]->neighbors.push_back(ps[j]);
          ps[j]->neighbors.push_back(ps[i]);
        }
      }
    }
    double ds = 0.0;
    cout << "avg rho: ";
    for (Particle *p : ps) {
      p->initializeWithNewNeighbors();
    }
    for (size_t i = 0; i < NEWTON_NUM_STEPS; i++) {
      for (Particle *p : ps) {
        p->newtonStepCalculateLambda();
      }
      if (i == 0) {
        for (Particle *p : ps) {
          ds += p->getLatestDensityEstimate();
        }
        cout << ds / ps.size() << " => ";
      }
      for (Particle *p : ps) {
        p->newtonStepUpdatePosition(bvh);
      }
    }
    for (Particle *p : ps) {
      p->updateVelocity(delta_t);
      p->calculateVorticityApplyXSPHViscosity();
    }
    ds = 0.0;
    for (Particle *p : ps) {
      p->applyVorticity(delta_t);
      ds += p->getLatestDensityEstimate();
      p->updatePosition();
    }
    cout << ds / ps.size() << endl;
    surfaceUpToTimestep = false;
  }

  void Particles::timeStep() {
    timeStep(DEFAULT_DELTA_T);
  }

  void Particles::redraw(const Color& c) {
    for (Particle *p : ps) {
      Misc::draw_sphere_opengl(p->getPosition(), VIS_RADIUS, p->getDensityBasedColor());
    }
  }

  GRIDCELL Particles::generate_gridcell(double x1, double x2, double y1, double y2, double z1, double z2) {
    GRIDCELL g;
    g.p[0] = Vector3D(x1, y1, z1);
    g.p[1] = Vector3D(x1, y2, z1);
    g.p[2] = Vector3D(x2, y2, z1);
    g.p[3] = Vector3D(x2, y1, z1);
    g.p[4] = Vector3D(x1, y1, z2);
    g.p[5] = Vector3D(x1, y2, z2);
    g.p[6] = Vector3D(x2, y2, z2);
    g.p[7] = Vector3D(x2, y1, z2);
    for (int i = 0; i < 8; i++) {
      g.val[i] = estimateDensityAt(g.p[i]);
    }
    return g;
  }

  /*
    hardcoded min and max coordinates of particels
    x: [-1,1]
    y: [0,1.5]
    z: [-1,1]
  */

  double Particles::get_min(int axis){
    switch(axis) {
      case(0):return -1;
      case(1):return 0;
      case(2):return -1;
      default: exit(-1);
    }
  }

  double Particles::get_max(int axis){
    switch(axis) {
      case(0):return 1;
      case(1):return 1.5;
      case(2):return 1;
      default: exit(-1);
    }
  }

  // void Particles::add_to_mesh(const Mesh* mesh, double fStepSize) {
  std::vector<Primitive *> Particles::getSurfacePrims(double isolevel, double fStepSize, BSDF* bsdf) {
    //, int &vertices_base, int &poly_base) {
    std::vector<Primitive *> prims;
    double xmin = get_min(0);
    double ymin = get_min(1);
    double zmin = get_min(2);
    double xmax = get_max(0);
    double ymax = get_max(1);
    double zmax = get_max(2);
    int xsteps = (int)((xmax-xmin)/fStepSize);
    int ysteps = (int)((ymax-ymin)/fStepSize);
    int zsteps = (int)((zmax-zmin)/fStepSize);

    vector<vector<Index> > polygons;
    vector<Vector3D> vertexPositions;

    for (int ix = 0; ix <= xsteps; ix++)
      for (int iy = 0; iy <= ysteps; iy++)
        for (int iz = 0; iz <= zsteps; iz++) {
          double x1 = xmin+ix*fStepSize;
          double x2 = x1+fStepSize;
          double y1 = ymin+iy*fStepSize;
          double y2 = y1+fStepSize;
          double z1 = zmin+iz*fStepSize;
          double z2 = z1+fStepSize;

          GRIDCELL grid = generate_gridcell(x1,x2,y1,y2,z1,z2);
          std::vector<TRIANGLE *> triangles = polygonise(grid, isolevel);
          for (int i=0; i<triangles.size(); i++) {
            Vector3D p1 = triangles[i]->p[0];
            Vector3D p2 = triangles[i]->p[1];
            Vector3D p3 = triangles[i]->p[2];
            Vector3D n1 = getVertexNormal(p1);
            Vector3D n2 = getVertexNormal(p2);
            Vector3D n3 = getVertexNormal(p3);

            MarchingTriangle *tri = new MarchingTriangle(p1, p2, p3, n1, n2, n3, bsdf);
            // cout << "p1:" << p1 << " p2:" << p2 << "p3:" << p3 << endl;
            // cout << "n1:" << n1 << " n2:" << n2 << "n3:" << n3 << endl;
            prims.push_back(tri);
            // add_triangle_to_mesh(mesh,p0,p1,p2);
          }
    }
    return prims;
  }

  void Particles::updateSurface() {
    if (surfaceUpToTimestep) {
      return;
    }
    surface = getSurfacePrims(ISO_LEVEL, H * 0.5,
      new DiffuseBSDF(Spectrum(0.1,0.1,0.8)));
    surfaceUpToTimestep = true;
  }


  //vGetNormal() finds the gradient of the scalar field at a point
  //This gradient can be used as a very accurate vertx normal for lighting calculations
  Vector3D Particles::getVertexNormal(Vector3D &pos) {
    double fX = pos[0];
    double fY = pos[1];
    double fZ = pos[2];
    Vector3D n;
    n.x = estimateDensityAt(Vector3D(fX-GRADIENT_EPS, fY, fZ)) - estimateDensityAt(Vector3D(fX+GRADIENT_EPS, fY, fZ));
    n.y = estimateDensityAt(Vector3D(fX, fY-GRADIENT_EPS, fZ)) - estimateDensityAt(Vector3D(fX, fY+GRADIENT_EPS, fZ));
    n.z = estimateDensityAt(Vector3D(fX, fY, fZ-GRADIENT_EPS)) - estimateDensityAt(Vector3D(fX, fY, fZ+GRADIENT_EPS));
    if (n.norm() > 0) {
      // cout << "n:" << n << endl;
      return n.unit();
    }
    return n;
  }


  string Particles::paramsString() {
    stringstream ss;
    ss << "Fluid simulation parameters: " << endl
      << "\tH: " << H << endl
      << "\tNewton steps: " << NEWTON_NUM_STEPS << endl
      << "\tConstraint relaxation epsilon: " << EPSILON << endl
      << "\tTensile artificial pressure coefficient K: " << K << endl
      << "\tTensile artificial pressure exponent N: " << N << endl
      << "\tVorticity confinement coefficient epsilon: " << VORTICITY_EPSILON << endl
      << "\tViscosity coefficient C: " << C << endl
      << "\tParticle visualization radius: " << VIS_RADIUS << endl
      << "\tCollision response velocity factor: " << VELOCITY_BOUNCE_FACTOR << endl;
    return ss.str();
  }

  void Particles::estimateDensities() {
    for (Particle *p : ps) {
      p->estimateDensityWithNeighbors(ps);
    }
  }

  double Particles::estimateDensityAt(Vector3D pos) {
    double density = 0.0;
    for (Particle *p : ps) {
        density += poly6_kernel(p->getPosition() - pos);
    }
    return density;

  }

}  // namespace CGL

#undef VIS_RADIUS
#undef ISO_LEVEL
#undef GRADIENT_EPS
#undef FSTEPSIZE
#undef VELOCITY_BOUNCE_FACTOR
#undef DEFAULT_DELTA_T
#undef H
#undef H2
// #undef MAX_NUM_NEIGHBORS
#undef NUM_NEIGHBOR_ALERT_THRESHOLD
#undef NEWTON_NUM_STEPS
#undef EPSILON
#undef K
#undef N
#undef VORTICITY_EPSILON
#undef C
// #undef CELL_NUM_X
// #undef CELL_NUM_Y
// #undef CELL_NUM_Z


