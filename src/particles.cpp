#include "particles.h"

using namespace CGL::StaticScene;

namespace CGL {

  void Particle::findNeighbors() {
    neighbors.clear();
    for (Particle *p : particles->ps) {
      if ((p->origin() - origin()).norm2() <= particles->sq_neighbor_radius) {
        neighbors.push_back(p);
      }
    }
    if (neighbors.size() < particles->p_num_neighbor_alert_thresh) {
      cerr << "Particle" << origin()
        << " only has " << neighbors.size()
        << " neighbors." << endl;
    }
  }

  void Particle::applyForces(double delta_t) {
    for (Force *f : particles->fs) {
      velocity += delta_t * f->getAccerlation(particles->simulate_time, *this);
    }
  }

  void Particle::applyVelocity(double delta_t) {
    origin() += delta_t * velocity;
  }

  bool Particles::simulateToTime(double t) {
    if (t <= simulate_time) {
      return false;
    }
    while (simulate_time + default_delta_t < t) {
      timeStep(default_delta_t);
    }
    if (t > simulate_time) {
      timeStep(t - simulate_time);
    }
    return true;
  }

  void Particles::timeStep(double delta_t) {
    simulate_time += delta_t;
  }

}  // namespace CGL
