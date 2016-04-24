#include "bbox.h"

#include "GL/glew.h"

#include <algorithm>
#include <iostream>

namespace CGL {

bool BBox::intersect(const Ray& r, double& t0, double& t1) const {

  // TODO Part 2, task 2:
  // Implement ray - bounding box intersection test
  // If the ray intersected the bounding box within the range given by
  // t0, t1, update t0 and t1 with the new intersection times.
  double tx0 = (min.x - r.o.x) / r.d.x;
  double tx1 = (max.x - r.o.x) / r.d.x;
  double ty0 = (min.y - r.o.y) / r.d.y;
  double ty1 = (max.y - r.o.y) / r.d.y;
  double tz0 = (min.z - r.o.z) / r.d.z;
  double tz1 = (max.z - r.o.z) / r.d.z;
  if (tx0 > tx1) {
    double t = tx0;
    tx0 = tx1;
    tx1 = t;
  }
  if (ty0 > ty1) {
    double t = ty0;
    ty0 = ty1;
    ty1 = t;
  }
  if (tz0 > tz1) {
    double t = tz0;
    tz0 = tz1;
    tz1 = t;
  }
  double tmin, tmax;
  if ((tx0 >= ty0) && (tx0 >= tz0)) {
    tmin = tx0;
  } else if (ty0 >= tz0) {
    tmin = ty0;
  } else {
    tmin = tz0;
  }
  if (tmin < t0) {
    tmin = t0;
  }
  if ((tx1 <= ty1) && (tx1 <= tz1)) {
    tmax = tx1;
  } else if (ty1 <= tz1) {
    tmax = ty1;
  } else {
    tmax = tz1;
  }
  if (tmax > t1) {
    tmax = t1;
  }
  if (tmax >= tmin) {
    t0 = tmin;
    t1 = tmax;
    return true;
  }
  return false;
}

void BBox::draw(Color c) const {

  glColor4f(c.r, c.g, c.b, c.a);

	// top
	glBegin(GL_LINE_STRIP);
	glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, max.y, min.z);
  glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, max.y, max.z);
  glVertex3d(max.x, max.y, max.z);
	glEnd();

	// bottom
	glBegin(GL_LINE_STRIP);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(min.x, min.y, max.z);
  glVertex3d(max.x, min.y, max.z);
  glVertex3d(max.x, min.y, min.z);
  glVertex3d(min.x, min.y, min.z);
	glEnd();

	// side
	glBegin(GL_LINES);
	glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, min.y, max.z);
	glVertex3d(max.x, max.y, min.z);
  glVertex3d(max.x, min.y, min.z);
	glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, min.y, min.z);
	glVertex3d(min.x, max.y, max.z);
  glVertex3d(min.x, min.y, max.z);
	glEnd();

}

std::ostream& operator<<(std::ostream& os, const BBox& b) {
  return os << "BBOX(" << b.min << ", " << b.max << ")";
}

} // namespace CGL
