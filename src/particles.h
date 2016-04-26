#ifndef CGL_PARTICLES_H
#define CGL_PARTICLES_H

#include <utility>
#include <algorithm>

#include "CGL/CGL.h"
#include "static_scene/sphere.h"

namespace CGL {

inline double poly6_kernel(Vector3D r, double h) {
  return (315 / (64 * PI * intpow<9>(h))) * intpow<3>(intpow<2>(h) - r.norm2());
}

inline Vector3D grad_spiky_kernel(Vector3D r, double h) {
  return (45 / (PI * intpow<6>(h))) * intpow<2>(h - r.norm()) * r.unit();
}

/* Allow use in CGL::StaticScene::Particle. */
struct Particles;

/* Keep the subclass in CGL::StaticScene as a good practice. */
namespace StaticScene {
  class Particle : public StaticScene::Sphere {
   public:
    Vector3D velocity;
    const double rest_density;
    const double mass;

    Particle(
      const StaticScene::SphereObject* object,
      const Vector3D& v,
      const Vector3D& pos,
      const double r,
      const double rest_density,
      const Particles *particles
    ) :
      StaticScene::Sphere(object, pos, r),
      velocity(v),
      rest_density(rest_density),
      particles(particles),
      mass(4.0 / 3 * PI * intpow<3>(r) * rest_density) { }

    bool collide(Particle *other);

    void findNeighbors();

    void applyForces(double delta_t);
    void applyVelocity(double delta_t);
    void calculateLambda();
    void esitmateDensity();

   private:
    const Particles *particles;
    double lambda; // changing during simulation
    double density; // changing during simulation
    std::vector<Particle *> neighbors; // changing during simulation
  };

} // namespace StaticScene

/* The above defined class. */
using CGL::StaticScene::Particle;

struct Force {
  virtual Vector3D getAccerlation(double t, Particle &p);
};

struct Gravity : Force {
  Vector3D getAccerlation(double t, Particle &p) {
    return Vector3D(0, 0, -9.8);
  }
};

struct Particles {
  std::vector<Particle *> ps;
  std::vector<Force *> fs;
  double simulate_time;

  // params
  double default_delta_t = 0.001;
  double sq_neighbor_radius = 4;
  // should be 30~40 for SPH density est to be stable
  double p_num_neighbor_alert_thresh = 18;
  // used in calculating lambda
  double epsilon = 0.1;

  Particles() : simulate_time(0.0) {};

  // return True iff t > current time.
  bool simulateToTime(double t);
  void timeStep(double delta_t);
};

} // namespace CGL

#endif // CGL_PARTICLES_H
