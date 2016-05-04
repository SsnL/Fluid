#include <sstream>

#include "particles.h"
#include "misc/sphere_drawing.h"

// params
// particle visualization radius
#define VIS_RADIUS 0.02
// the threshold of isovalue for isosurface
#define ISO_LEVEL 0.8
// the resolution of isosurface
#define FSTEPSIZE 0.25
// default time step
#define DEFAULT_DELTA_T 0.016
// used in SPH density estimation
#define H 0.5
// h^2, neighbor radius squared
#define H2 0.25
// max num neighbors
// #define MAX_NUM_NEIGHBORS 45
// should be 30~40 for SPH density est to be stable
#define NUM_NEIGHBOR_ALERT_THRESHOLD 18
// number of Newton steps
#define NEWTON_NUM_STEPS 12
// used in calculating lambda
#define EPSILON 5
// artificial pressure coeff
#define K 0.001
// artificial pressure exp
#define N 4
// vorticity
#define VORTICITY_EPSILON 0.01
// XSPH viscocity
#define C 0.0001

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
    // to avoid particles stuck between light and ceiling, clip it
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
    velocity.y -= 9.8 * delta_t;
    new_position = position;
    clamp(new_position, velocity * delta_t, bvh);
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
    velocity += delta_t * VORTICITY_EPSILON * cross(grad_vorticity.unit(), vorticity);
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
    meshUpToTimestep = false;
  }

  void Particles::timeStep() {
    timeStep(DEFAULT_DELTA_T);
  }

  void Particles::redraw(const Color& c) {
    for (Particle *p : ps) {
      Misc::draw_sphere_opengl(p->getPosition(), VIS_RADIUS, p->getDensityBasedColor());
    }
  }

  double Particles::estimateDensityAt(Vector3D pos) {
    return 0;
  }

  GRIDCELL Particles::generate_gridcell(float x1, float x2, float y1, float y2, float z1, float z2) {
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


  void Particles::add_triangle(const vector<vector<Index> >& polygons,
                               const vector<Vector3D>& vertexPositions,
                               Vector3D p0, Vector3D p1, Vector3D p2) {

     size_t base = mesh->vertices.size();
     mesh->vertices.push_back(v0);
     mesh->vertices.push_back(v1);
     mesh->vertices.push_back(v2);
     Collada::Polygon poly;
     poly.vertex_indices.push_back(base);
     poly.vertex_indices.push_back(base + 1);
     poly.vertex_indices.push_back(base + 2);
     mesh->polygons.push_back(poly);
  }

  // void Particles::remove_from_mesh(const Mesh* mesh, int vertices_base, int poly_base) {
  //    int vertices_pops = mesh->vertices.size()-vertices_base;
  //    int poly_pops = mesh->polygons.size()-vertices_base;
  //    for (int i = 0; i < vertices_pops; i++) {
  //     mesh->vertices.pop_back();
  //    }
  //    for (int i = 0; i < poly_pops; i++) {
  //     mesh->polygons.pop_back();
  //    }
  // }

  // void Particles::add_to_mesh(const Mesh* mesh, float fStepSize) {
  DynamicScene::Mesh* Particles::get_mesh(float fStepSize) {
    //, int &vertices_base, int &poly_base) {
    DynamicScene::Mesh *mesh = new DynamicScene::Mesh();
    double isolevel; // the threshold isolevel
    float xmin = get_min(0);
    float ymin = get_min(1);
    float zmin = get_min(2);
    float xmax = get_max(0);
    float ymax = get_max(1);
    float zmax = get_max(2);
    int xsteps = (int)((xmax-xmin)/fStepSize);
    int ysteps = (int)((ymax-ymin)/fStepSize);
    int zsteps = (int)((zmax-zmin)/fStepSize);
    // vertices_base = mesh->vertices.size();
    // poly_base = mesh->polygons.size();
    
    vector<vector<Index> > polygons;
    vector<Vector3D> vertexPositions;
    
    for (int ix = 0; ix <= xsteps; ix++) 
      for (int iy = 0; iy <= ysteps; iy++) 
        for (int iz = 0; iz <= zsteps; iz++) {
          float x1 = xmin+ix*fStepSize;
          float x2 = x1+fStepSize;
          float y1 = ymin+iy*fStepSize;
          float y2 = y1+fStepSize;
          float z1 = zmin+iz*fStepSize;
          float z2 = z1+fStepSize;
          GRIDCELL grid = generate_gridcell(x1,x2,y1,y2,z1,z2);
          TRIANGLE *triangles = TRIANGLE[16]; 
          int ntriangles = Polygonise(grid, isolevel, triangles);
          for (int i=0; i<ntriangles; i++) {
            Vector3D p0 = triangles[i][0];
            Vector3D p1 = triangles[i][1];
            Vector3D p2 = triangles[i][2];
            add_triangle_to_mesh(mesh,p0,p1,p2);
          }
    }
    return mesh;
  }

  DynamicScene::Mesh *updateMesh() {
    if (meshUpToTimestep) {
      return mesh;
    }
    // int vertices_base, poly_base;
    // mesh = new Mesh();
    // add_to_mesh(mesh, FSTEPSIZE, vertices_base, poly_base);
    mesh = get_mesh(FSTEPSIZE);
    meshUpToTimestep = true;
    return mesh;
  }


  //vGetNormal() finds the gradient of the scalar field at a point
  //This gradient can be used as a very accurate vertx normal for lighting calculations
  Vector3D Particles::vGetNormal(Vector3D &pos)
  {
    double fX = pos[0];
    double fY = pos[1];
    double fZ = pos[2];
    Vector3D n;
    n.x = estimateDensityAt(Vector3D(fX-EPS_D, fY, fZ)) - estimateDensityAt(Vector3D(fX+EPS_D, fY, fZ));
    n.y = estimateDensityAt(Vector3D(fX, fY-EPS_D, fZ)) - estimateDensityAt(Vector3D(fX, fY+EPS_D, fZ));
    n.z = estimateDensityAt(Vector3D(fX, fY, fZ-EPS_D)) - estimateDensityAt(Vector3D(fX, fY, fZ+EPS_D));
    return n.unit();
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
      << "\tParticle visualization radius: " << VIS_RADIUS << endl;
    return ss.str();
  }

}  // namespace CGL

#undef VIS_RADIUS
#undef ISO_LEVEL
#undef FSTEPSIZE

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


