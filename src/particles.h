#ifndef CGL_PARTICLES_H
#define CGL_PARTICLES_H

#include <utility>
#include <algorithm>

#include "CGL/CGL.h"
#include "static_scene/sphere.h"
#include "bvh.h"

// params
// used in SPH density estimation
#define H 2

using namespace CGL::StaticScene;

namespace CGL {

inline double poly6_kernel(Vector3D r) {
  double r_l = r.norm();
  if (r_l > H) {
    return 0.0;
  }
  double temp = intpow<2>(H) - intpow<2>(r_l);
  return (315 * intpow<3>(temp)) / (64 * PI * intpow<9>(H));
}

inline Vector3D grad_spiky_kernel(Vector3D r) {
  double r_l = r.norm();
  if (r_l > H || r_l == 0.0) {
    return Vector3D();
  }
  return (45 * intpow<2>(H - r_l) * r) / (PI * intpow<6>(H) * r_l);
}

/* Allow use in CGL::StaticScene::Particle. */
struct Particles;

/* Keep the subclass in CGL::StaticScene as a good practice. */
namespace StaticScene {
  class Particle : public Sphere {
   public:
    Vector3D velocity;
    const double rest_density;
    // const double mass;
    Vector3D new_origin;
    std::vector<Particle *> neighbors; // changing during simulation
    // Vector3D vorticity;

    static double tensile_instability_scale;

    Particle(
      const SphereObject* object,
      const Vector3D& v,
      const double rest_density
    ) :
      StaticScene::Sphere(object, object->o, object->r),
      velocity(v),
      rest_density(rest_density) { }

    double getLatestDensityEstimate() {
      return density;
    }

    void initializeWithNewNeighbors();

    // naive force and velocity
    void applyForceVelocity(double delta_t);

    // Newton step density constraint
    void newtonStepCalculateLambda();
    void newtonStepUpdatePosition(BVHAccel *bvh);

    // update velocity basing on new position
    // positions is unchanged
    void updateVelocity(double delta_t);

    // calculate vorticity and apply XSPH viscosity
    void calculateVorticityApplyXSPHViscosity();

    // apply vorticity confinement
    void applyVorticity(double delta_t);

    // update position to
    void updatePosition();

   private:
    // const Particles *particles;
    // SPH estimate, changing during simulation
    double density;
    double lambda; // changing during simulation
    Vector3D vorticity; // changing during simulation
    std::vector<double> s_corr; // changing during simulation, same order as neighbors
    std::vector<Vector3D> grad_w_neighbors; // changing during simulation, same order as neighbors
  };

} // namespace StaticScene

/* The above defined class. */
// using CGL::StaticScene::Particle;

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
  // std::vector<Force *> fs;
  BVHAccel* bvh;
  double simulate_time;

  // TODO: Also set PS, FS, and BVH
  Particles() : bvh(NULL), simulate_time(0.0) {
    for (int i = -2; i < 3; i++)
      for (int j = 0; j < 3; j++)
        for (int k = -2; k < 3; k++)
          ps.push_back(new Particle(
            new StaticScene::SphereObject(Vector3D(0.12 * i, 0.12 * j + 1, 0.12 * k), 0.05f, new DiffuseBSDF(Spectrum(0.5f,0.5f,0.5f))),
            Vector3D(0, 0, 2),
            1.0f
          ));
  };

  // return True iff t > current time.
  bool simulateToTime(double t);
  void timeStep(double delta_t);
  void redraw(const Color& c);
};

} // namespace CGL

#undef H

#endif // CGL_PARTICLES_H
