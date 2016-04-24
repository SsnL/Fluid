#include "sphere.h"

#include <cmath>

#include  "../bsdf.h"
#include "../misc/sphere_drawing.h"

namespace CGL { namespace StaticScene {

bool Sphere::test(const Ray& r, double& t1, double& t2) const {

  // TODO Part 1, task 4:
  // Implement ray - sphere intersection test.
  // Return true if there are intersections and writing the
  // smaller of the two intersection times in t1 and the larger in t2.
  double a = dot(r.d, r.d);
  Vector3D m = r.o - o;
  double b = 2 * dot(m, r.d);
  double c = dot(m, m) - r2;
  double delta = b * b - 4 * a * c;
  if (delta < 0) {
    return false;
  }
  double t = (-b - sqrt(delta)) / (2 * a);
  if ((t >= r.min_t) && (t <= r.max_t)) {
    t1 = t;
    t = (-b + sqrt(delta)) / (2 * a);
    if ((t >= r.min_t) && (t <= r.max_t)) {
      t2 = t;
    }
  } else {
    t = (-b + sqrt(delta)) / (2 * a);
    if ((t >= r.min_t) && (t <= r.max_t)) {
      t1 = t;
    } else {
      return false;
    }
  }
  return true;

}

bool Sphere::intersect(const Ray& r) const {

  // TODO Part 1, task 4:
  // Implement ray - sphere intersection.
  // Note that you might want to use the the Sphere::test helper here.
  double t1, t2;
  if (test(r, t1, t2)) {
    r.max_t = t1;
    return true;
  }
  return false;

}

bool Sphere::intersect(const Ray& r, Intersection *i) const {

  // TODO Part 1m task 4:
  // Implement ray - sphere intersection.
  // Note again that you might want to use the the Sphere::test helper here.
  // When an intersection takes place, the Intersection data should be updated
  // correspondingly.
  double t1, t2;
  if (test(r, t1, t2)) {
    r.max_t = t1;
    i->t = t1;
    Vector3D n = r.o + t1 * r.d - o;
    i->n = n.unit();
    i->primitive = this;
    i->bsdf = get_bsdf();
    return true;
  }
  return false;

}

void Sphere::draw(const Color& c) const {
  Misc::draw_sphere_opengl(o, r, c);
}

void Sphere::drawOutline(const Color& c) const {
    //Misc::draw_sphere_opengl(o, r, c);
}


} // namespace StaticScene
} // namespace CGL
