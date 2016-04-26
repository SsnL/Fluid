#ifndef CGL_PARTICLES_H
#define CGL_PARTICLES_H

#include <utility>
#include <algorithm>

#include "CGL/CGL.h"
#include "static_scene/sphere.h"

namespace CGL {

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
      const double rd,
      const Particles *particles
    ) :
      StaticScene::Sphere(object, pos, r),
      velocity(v),
      rest_density(rd),
      particles(particles),
      mass(4.0 / 3 * PI * intpow<3>(r) * rest_density) { }

    bool collide(Particle *other);

   private:
    const Particles *particles;
    double lambda;
    std::vector<Particle *> neighbors;
    void findNeighbors();
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
  double default_delta_t = 0.02;
  double sq_neighbor_radius = 4;

  Particles() : simulate_time(0.0) {};

  // return True iff t > current time.
  bool simulateToTime(double t);
  void timeStep(double delta_t);
};

} // namespace CGL

#endif // CGL_PARTICLES_H
