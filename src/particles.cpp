#include "particles.h"
#include "misc/sphere_drawing.h"
#include "GL/glew.h"

// params
// artificial pressure denom
#define DEFAULT_DELTA_T 0.016
// used in SPH density estimation
#define H 2
// h^2, neighbor radius squared
#define H2 4
// should be 30~40 for SPH density est to be stable
#define NUM_NEIGHBOR_ALERT_THRESHOLD 18
// number of Newton steps
#define NEWTON_NUM_STEPS 7
// used in calculating lambda
#define EPSILON 0.1
// artificial pressure coeff
#define K 0.1
// artificial pressure exp
#define N 4
// vorticity
#define VORTICITY_EPSILON 5
// XSPH viscocity
#define C 0.01

using namespace CGL::StaticScene;

namespace CGL {
  double Particle::tensile_instability_scale = 1.0 / poly6_kernel(Vector3D(0, 0, 0.2 * H));

  void Particle::initializeWithNewNeighbors() {
    s_corr.resize(neighbors.size());
    grad_w_neighbors.resize(neighbors.size());
    if (neighbors.size() < NUM_NEIGHBOR_ALERT_THRESHOLD) {
      cerr << "Particle" << origin()
        << " only has " << neighbors.size()
        << " neighbors." << endl;
    }
  }

  void Particle::applyForceVelocity(double delta_t) {
    // for (Force *f : particles->fs) {
    //   velocity += delta_t * f->getAccerlation(particles->simulate_time, *this);
    // }
    // assume only gravity for now
    velocity.z -= 9.8 * delta_t;
    new_origin = origin() + delta_t * velocity;
  }

  void Particle::newtonStepCalculateLambda() {
    density = 0.0;
    double denom = 0.0;
    Vector3D grad_i_c_i;
    for (size_t i = 0; i < neighbors.size(); i++) {
      Vector3D r = new_origin - neighbors[i]->new_origin;
      double w_i_n = poly6_kernel(r);
      Vector3D grad_w_i_n = grad_spiky_kernel(r);
      // SPH density esitmator
      density += w_i_n;
      grad_i_c_i += grad_w_i_n;
      denom += grad_w_i_n.norm2();
      s_corr[i] = -K * intpow<N>(w_i_n / tensile_instability_scale);
      grad_w_neighbors[i] = grad_w_i_n;
    }
    double c_i = density / rest_density - 1.0;
    denom += grad_i_c_i.norm2();
    lambda = -c_i / (denom / intpow<2>(rest_density) + EPSILON);
  }

  void Particle::newtonStepUpdatePosition(BVHAccel *bvh) {
    Vector3D delta_p;
    for (size_t i = 0; i < neighbors.size(); i++) {
      delta_p += (neighbors[i]->lambda + lambda + s_corr[i]) * grad_w_neighbors[i];
    }
    delta_p /= rest_density;
    Vector3D direction = delta_p.unit();
    Ray r(origin(), direction, delta_p.norm());
    bvh->intersect(r);
    delta_p = direction * r.max_t;
    new_origin += delta_p;
  }

  void Particle::updateVelocity(double delta_t) {
    velocity = (new_origin - origin()) / delta_t;
  }

  void Particle::calculateVorticityApplyXSPHViscosity() {
    vorticity = Vector3D();
    Vector3D viscocity;
    for (Particle *n : neighbors) {
      Vector3D v_i_n = n->velocity - velocity;
      Vector3D r = new_origin - n->new_origin;
      double w_i_n = poly6_kernel(r);
      Vector3D grad_w_i_n = grad_spiky_kernel(r);
      vorticity += cross(grad_w_i_n, v_i_n);
      viscocity += v_i_n * w_i_n;
    }
    velocity += C * viscocity;
  }

  void Particle::applyVorticity(double delta_t) {
    Vector3D grad_vorticity;
    for (Particle *n : neighbors) {
      Vector3D r = new_origin - n->new_origin;
      grad_vorticity += n->vorticity * poly6_kernel(r);
    }
    velocity += delta_t * VORTICITY_EPSILON * cross(grad_vorticity.norm(), vorticity);
  }

  void Particle::updatePosition() {
    origin() = new_origin;
  }

  bool Particles::simulateToTime(double t) {
    if (t <= simulate_time) {
      return false;
    }
    while (simulate_time + DEFAULT_DELTA_T < t) {
      timeStep(DEFAULT_DELTA_T);
    }
    if (t > simulate_time) {
      timeStep(t - simulate_time);
    }
    return true;
  }

  void Particles::timeStep(double delta_t) {
    simulate_time += delta_t;
    for (size_t i = 0; i < ps.size(); i++) {
      for (size_t j = i + 1; j < ps.size(); j++) {
        if ((ps[i]->origin() - ps[j]->origin()).norm2() <= H2) {
          ps[i]->neighbors.push_back(ps[j]);
          ps[j]->neighbors.push_back(ps[i]);
        }
      }
    }
    for (Particle *p : ps) {
      p->initializeWithNewNeighbors();
      p->applyForceVelocity(delta_t);
    }
    for (size_t i = 0; i < NEWTON_NUM_STEPS; i++) {
      for (Particle *p : ps) {
        p->newtonStepCalculateLambda();
      }
      for (Particle *p : ps) {
        p->newtonStepUpdatePosition(bvh);
      }
    }
    for (Particle *p : ps) {
      p->updateVelocity(delta_t);
      p->calculateVorticityApplyXSPHViscosity();
    }
    for (Particle *p : ps) {
      p->applyVorticity(delta_t);
      p->updatePosition();
    }
  }

  void Particles::redraw(const Color& c) {
    // simulateToTime(simulate_time + DEFAULT_DELTA_T);
    // glPolygonOffset(1.0, 1.0);
    // glEnable(GL_POLYGON_OFFSET_FILL);
    // Misc::draw_sphere_opengl(Vector3D(), 1, Color(1,1,1,1));
    for (Particle *p : ps) {
      Misc::draw_sphere_opengl(p->origin(), p->radius(), c);
    }
    // glDisable(GL_POLYGON_OFFSET_FILL);
  }

}  // namespace CGL

#undef DEFAULT_DELTA_T
#undef H
#undef H2
#undef NUM_NEIGHBOR_ALERT_THRESHOLD
#undef NEWTON_NUM_STEPS
#undef EPSILON
#undef K
#undef N
#undef VORTICITY_EPSILON
#undef C


