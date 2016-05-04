#ifndef CGL_PARTICLES_H
#define CGL_PARTICLES_H

#include <utility>
#include <algorithm>

#include "CGL/CGL.h"
#include "bsdf.h"
#include "static_scene/sphere.h"
#include "bvh.h"
#include "random_util.h"
#include "dynamic_scene/mesh.h"
#include "marching.h"

using namespace CGL::StaticScene;

namespace CGL {

class Particle {
 public:
  Vector3D velocity;
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

  inline double const getLatestDensityEstimate() {
    return density;
  }

  inline const Vector3D getPosition() {
    return position;
  }

  inline const Vector3D getNewPosition() {
    return new_position;
  }

  Color getDensityBasedColor() {
    double ratio = (density - 0.8 * rest_density) / (0.4 * rest_density),
      clamp_ratio = min(1.0, max(0.0, ratio));
    return Color(1.0, 1.0 - clamp_ratio, 1.0 - clamp_ratio, 1.0);
  }

  // given a set of neighbors, estimate density
  void estimateDensityWithNeighbors(std::vector<Particle *> neighbors);

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
  const double rest_density;
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
  const double rest_density;
  std::vector<Primitive *> surface;
  bool surfaceUpToTimestep = false;

  Particles(double rest_density = 1000.0) : bvh(NULL), rest_density(rest_density), simulate_time(0.0) {
    estimateDensities();
  };

  void addParticle(Vector3D pos, Vector3D v) {
    ps.push_back(new Particle(pos, v, rest_density));
  }

  // fluid simulation
  void timeStep(double delta_t);
  void timeStep();
  void estimateDensities();
  double estimateDensityAt(Vector3D pos);

  // marching cube
  void updateSurface();
  Vector3D getVertexNormal(Vector3D &pos);
  GridCell generateGridCell(double x1, double x2, double y1, double y2, double z1, double z2);
  double getParticlesAxisMin(int axis);
  double getParticlesAxisMax(int axis);
  std::vector<Primitive *> getSurfacePrims(double isolevel, double fStepSize, BSDF* bsdf);

  // utils
  string paramsString();
  void redraw(const Color& c);

};

} // namespace CGL

#endif // CGL_PARTICLES_H
