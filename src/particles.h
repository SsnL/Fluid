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

struct Particles;

class Particle {
 public:
  Vector3D velocity;
  const double rest_density;
  std::vector<Particle *> neighbors; // changing during simulation

  static double tensile_instability_scale;

  Particle(
    const Vector3D &p,
    const Vector3D &v,
    const double rest_density
  ) :
    position(p),
    new_position(p),
    velocity(v),
    rest_density(rest_density) { }

  double getLatestDensityEstimate() {
    return density;
  }

  inline Vector3D getPosition() {
    return position;
  }

  inline Vector3D getNewPosition() {
    return new_position;
  }

  Color getDensityBasedColor() {
    double ratio = (density - 0.8 * rest_density) / (0.4 * rest_density),
      clamp_ratio = min(1.0, max(0.0, ratio));
    return Color(1.0, 1.0 - clamp_ratio, 1.0 - clamp_ratio, 1.0);
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
  // accurate between simulation timesteps
  Vector3D position;
  // always ``newer'' than positions in simulation process sense
  Vector3D new_position;
  // SPH estimate, changing during simulation
  double density;
  double lambda; // changing during simulation
  Vector3D vorticity; // changing during simulation
  std::vector<double> s_corr; // changing during simulation, same order as neighbors
  std::vector<Vector3D> grad_w_neighbors; // changing during simulation, same order as neighbors
};

std::ostream &operator<<(std::ostream &os, const Particle &v);

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

  Particles() : bvh(NULL), simulate_time(0.0) {
    // for (int i = -6; i < 6; i++)
    //   for (int j = 0; j < 5; j++)
    //     for (int k = -6; k < 6; k++)
    //       ps.push_back(new Particle(
    //         Vector3D(0.15 * i, 0.2 * j + 0.5, 0.15 * k),
    //         Vector3D(0, -0.01, 0),
    //         150.0f
    //       ));
  };

  void timeStep(double delta_t);
  void timeStep();
  void redraw(const Color& c);
  string paramsString();
};

} // namespace CGL

#endif // CGL_PARTICLES_H
