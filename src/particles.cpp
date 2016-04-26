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
