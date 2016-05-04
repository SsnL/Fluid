#ifndef CGL_PARTICLES_H
#define CGL_PARTICLES_H

#include <utility>
#include <algorithm>

#include "CGL/CGL.h"
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
  DynamicScene::Mesh* mesh = new DynamicScene::Mesh();
  bool meshUpToTimestep = false;

  Particles() : bvh(NULL), simulate_time(0.0) {
    for (int i = 1; i < 10; i++)
      for (int j = 1; j < 14; j++)
        for (int k = 1; k < 10; k++)
          ps.push_back(new Particle(
            Vector3D(0.1 * i - 1, 0.1 * j, 1 - 0.1 * k),
            Vector3D(0, -1, 0),
            700.0f
          ));
    for (int i = 1; i < 10; i++)
      for (int j = 1; j < 14; j++)
        for (int k = 1; k < 10; k++)
          ps.push_back(new Particle(
            Vector3D(1 - 0.1 * i, 0.1 * j, 0.1 * k - 1),
            Vector3D(0, -1, 0),
            700.0f
          ));
  };

  void timeStep(double delta_t);
  void timeStep();
  void redraw(const Color& c);
  double estimateDensityAt(Vector3D pos);
  DynamicScene::Mesh *updateMesh();
  Vector3D vGetNormal(Vector3D &pos);
  string paramsString();

  GRIDCELL generate_gridcell(float x1, float x2, float y1, float y2, float z1, float z2);
  double get_min(int axis);
  double get_max(int axis);
  void add_triangle_to_mesh(const DynamicScene::Mesh* mesh, Vector3D p0, Vector3D p1, Vector3D p2);
  DynamicScene::Mesh* get_mesh(float fStepSize);
  
};

} // namespace CGL

#endif // CGL_PARTICLES_H
