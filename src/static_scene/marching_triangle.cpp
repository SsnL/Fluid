#include "marching_triangle.h"

#include "CGL/CGL.h"
#include "GL/glew.h"

namespace CGL { namespace StaticScene {

MarchingTriangle::MarchingTriangle(Vector3D p1, Vector3D p2, Vector3D p3,
  Vector3D n1, Vector3D n2, Vector3D n3, BSDF *bsdf) :
    p1(p2), p2(p3), p3(p3), n1(n1), n2(n2), n3(n3), bsdf(bsdf) { }

BBox MarchingTriangle::get_bbox() const {

  BBox bb(p1);
  bb.expand(p2);
  bb.expand(p3);
  return bb;

}

bool MarchingTriangle::intersect(const Ray& r) const {

  Vector3D e1 = p2 - p1;
  Vector3D e2 = p3 - p1;
  Vector3D s = r.o - p1;
  Vector3D s1 = cross(r.d, e2);
  Vector3D s2 = cross(s, e1);
  double d = dot(s1, e1);
  if (d == 0) {
    return false;
  }
  double t = dot(s2, e2) / d;
  if ((t < r.min_t) || (t > r.max_t)) {
    return false;
  }
  double u = dot(s1, s) / d;
  double v = dot(s2, r.d) / d;
  double w = 1 - u - v;
  if ((u < 0) || (u > 1) || (v < 0) || (v > 1) || (w < 0) || (w > 1)) {
    return false;
  }
  r.max_t = t;
  return true;
}

bool MarchingTriangle::intersect(const Ray& r, Intersection *isect) const {

  Vector3D e1 = p2 - p1;
  Vector3D e2 = p3 - p1;
  Vector3D s = r.o - p1;
  Vector3D s1 = cross(r.d, e2);
  Vector3D s2 = cross(s, e1);
  double d = dot(s1, e1);
  if (d == 0) {
    return false;
  }
  double t = dot(s2, e2) / d;
  if ((t < r.min_t) || (t > r.max_t)) {
    return false;
  }
  double u = dot(s1, s) / d;
  double v = dot(s2, r.d) / d;
  double w = 1 - u - v;
  if ((u < 0) || (u > 1) || (v < 0) || (v > 1) || (w < 0) || (w > 1)) {
    return false;
  }
  r.max_t = t;
  isect->t = t;
  isect->n = w * n1 + u * n2 + v * n3;
  isect->primitive = this;
  isect->bsdf = get_bsdf();
  return true;
}

void MarchingTriangle::draw(const Color& c) const {
  glColor4f(c.r, c.g, c.b, c.a);
  glBegin(GL_TRIANGLES);
  glVertex3d(p1.x, p1.y, p1.z);
  glVertex3d(p2.x, p2.y, p2.z);
  glVertex3d(p3.x, p3.y, p3.z);
  glEnd();
}

void MarchingTriangle::drawOutline(const Color& c) const {
  glColor4f(c.r, c.g, c.b, c.a);
  glBegin(GL_LINE_LOOP);
  glVertex3d(p1.x, p1.y, p1.z);
  glVertex3d(p2.x, p2.y, p2.z);
  glVertex3d(p3.x, p3.y, p3.z);
  glEnd();
}



} // namespace StaticScene
} // namespace CGL
