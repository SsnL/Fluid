#ifndef CGL_STATICSCENE_TRIANGLE_H
#define CGL_STATICSCENE_TRIANGLE_H

#include "object.h"
#include "primitive.h"

namespace CGL { namespace StaticScene {

/**
 * A single triangle from a mesh.
 * To save space, it holds a pointer back to the data in the original mesh
 * rather than holding the data itself. This means that its lifetime is tied
 * to that of the original mesh. The primitive may refer back to the mesh
 * object for other information such as normal, texcoord, material.
 */
 class MarchingTriangle : public Primitive {
 public:

  /**
   * Constructor.
   * \param p1 triangle vertex
   * \param p2 triangle vertex
   * \param p3 triangle vertex
   */
  MarchingTriangle(Vector3D p1, Vector3D p2, Vector3D p3,
    Vector3D n1, Vector3D n2, Vector3D n3, BSDF *bsdf);

   /**
    * Get the world space bounding box of the triangle.
    * \return world space bounding box of the triangle
    */
  BBox get_bbox() const;

   /**
    * Ray - Triangle intersection.
    * Check if the given ray intersects with the triangle, no intersection
    * information is stored.
    * \param r ray to test intersection with
    * \return true if the given ray intersects with the triangle,
              false otherwise
    */
  bool intersect(const Ray& r) const;

   /**
    * Ray - Triangle intersection 2.
    * Check if the given ray intersects with the triangle, if so, the input
    * intersection data is updated to contain intersection information for the
    * point of intersection.
    * \param r ray to test intersection with
    * \param i address to store intersection info
    * \return true if the given ray intersects with the triangle,
              false otherwise
    */
  bool intersect(const Ray& r, Intersection* i) const;

  /**
   * Get BSDF.
   */
  BSDF *get_bsdf() const { return bsdf; }

  /**
   * Draw with OpenGL (for visualizer)
   */
  void draw(const Color& c) const;

  /**
   * Draw outline with OpenGL (for visualizer)
   */
  void drawOutline(const Color& c) const;

 private:

  const Vector3D p1;
  const Vector3D p2;
  const Vector3D p3;
  const Vector3D n1;
  const Vector3D n2;
  const Vector3D n3;
  BSDF *bsdf;

}; // class Triangle

} // namespace StaticScene
} // namespace CGL

#endif //CGL_STATICSCENE_TRIANGLE_H
