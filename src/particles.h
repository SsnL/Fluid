#ifndef CGL_PARTICLES_H
#define CGL_PARTICLES_H

#include <utility>
#include <algorithm>

#include "CGL/CGL.h"
#include "static_scene/sphere.h"
#include "bvh.h"
#include "random_util.h"

using namespace CGL::StaticScene;

namespace CGL {

/* Allow use in CGL::StaticScene::Particle. */
struct Particles;

/* Keep the subclass in CGL::StaticScene as a good practice. */
namespace StaticScene {
  class Particle : public Sphere {
   public:
    Vector3D velocity;
    const double rest_density;
    // const double mass;
    std::vector<Particle *> neighbors; // changing during simulation
    // Vector3D vorticity;
    Color color;

    static double tensile_instability_scale;

    Particle(
      const SphereObject* object,
      const Vector3D& v,
      const double rest_density,
      Color c = Color(1, 1, 1, 1)
    ) :
      StaticScene::Sphere(object, object->o, object->r),
      velocity(v),
      rest_density(rest_density),
      color(c),
      new_origin(object->o) { }

    double getLatestDensityEstimate() {
      return density;
    }

    Vector3D newOrigin() {
      return new_origin;
    }

    void initializeWithNewNeighbors();

    // naive force and velocity
    void applyForceVelocity(BVHAccel *bvh, double delta_t);

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
    Vector3D new_origin;
    // SPH estimate, changing during simulation
    double density;
    double lambda; // changing during simulation
    Vector3D vorticity; // changing during simulation
    std::vector<double> s_corr; // changing during simulation, same order as neighbors
    std::vector<Vector3D> grad_w_neighbors; // changing during simulation, same order as neighbors
  };

} // namespace StaticScene

std::ostream& operator<<( std::ostream& os, const Particle& v );

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
    // for (int i = -6; i < 6; i++)
    //   for (int j = 0; j < 5; j++)
    //     for (int k = -6; k < 6; k++)
    //       ps.push_back(new Particle(
    //         new StaticScene::SphereObject(Vector3D(0.15 * i, 0.2 * j + 0.5, 0.15 * k), 0.02f, NULL),
    //         Vector3D(0, -0.01, 0),
    //         150.0f,
    //         Color(random_uniform(), random_uniform(), random_uniform(), 1)
    //       ));
  };

  void timeStep(double delta_t);
  void timeStep();
  void redraw(const Color& c);
};

} // namespace CGL

#endif // CGL_PARTICLES_H
