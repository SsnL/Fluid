#ifndef CGL_PARTICLES_H
#define CGL_PARTICLES_H

#include <utility>
#include <algorithm>

#include "CGL/CGL.h"
#include "static_scene/sphere.h"
#include "pathtracer.h"

namespace CGL {

class Particle : public StaticScene::Sphere {
 public:
  Vector3D velocity;
  double rest_density;
  double lambda;
  double mass;
  // std::vector<particle> neighbors;

  Particle(const StaticScene::SphereObject* object, const Vector3D& v, const Vector3D& pos, const double r,
    const double rd) : StaticScene::Sphere(object, pos, r), velocity(v), rest_density(rd) {
    mass = 4.0 / 3 * PI * r * r * r * rest_density;
  }

  bool collide(Particle *other);

};

struct Force {
  virtual Vector3D getAccerlation(double t, Particle &p);
};

struct Gravity : Force {
  Vector3D getAccerlation(double t, Particle &p) {
    return Vector3D(0, 0, -9.8 * p.mass);
  }
};

struct Particles {
  std::vector<Particle *> ps;
  std::vector<Force *> fs;
  PathTracer* pt;

  void simulateToTime(double t) {
    //
  };
};

}

// pt->time, pt->particles

#endif // CGL_PARTICLES_H
