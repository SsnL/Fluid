#include "bvh.h"

#include "CGL/CGL.h"
#include "static_scene/triangle.h"

#include <iostream>
#include <stack>

using namespace std;

namespace CGL { namespace StaticScene {

BVHAccel::BVHAccel(const std::vector<Primitive *> &_primitives,
                   size_t max_leaf_size) {

  root = construct_bvh(_primitives, max_leaf_size);

}

BVHAccel::~BVHAccel() {
  if (root) delete root;
}

BBox BVHAccel::get_bbox() const {
  return root->bb;
}

void BVHAccel::draw(BVHNode *node, const Color& c) const {
  if (node->isLeaf()) {
    for (Primitive *p : *(node->prims))
      p->draw(c);
  } else {
    draw(node->l, c);
    draw(node->r, c);
  }
}

void BVHAccel::drawOutline(BVHNode *node, const Color& c) const {
  if (node->isLeaf()) {
    for (Primitive *p : *(node->prims))
      p->drawOutline(c);
  } else {
    drawOutline(node->l, c);
    drawOutline(node->r, c);
  }
}

BVHNode *BVHAccel::construct_bvh(const std::vector<Primitive*>& prims, size_t max_leaf_size) {
  
  // TODO Part 2, task 1:
  // Construct a BVH from the given vector of primitives and maximum leaf
  // size configuration. The starter code build a BVH aggregate with a
  // single leaf node (which is also the root) that encloses all the
  // primitives.


  BBox centroid_box, bbox;

  for (Primitive *p : prims) {
    BBox bb = p->get_bbox();
    bbox.expand(bb);
    Vector3D c = bb.centroid();
    centroid_box.expand(c);
  }

  // You'll want to adjust this code.
  // Right now we just return a single node containing all primitives.
  BVHNode *node = new BVHNode(bbox);
  if (prims.size() <= max_leaf_size) {
    node->prims = new vector<Primitive *>(prims);
    return node;
  }
  Vector3D extent = bbox.extent;
  int split_dim;
  if ((extent.x >= extent.y) && (extent.x >= extent.z)) {
    split_dim = 0;
  } else if (extent.y >= extent.z) {
    split_dim = 1;
  } else {
    split_dim = 2;
  }
  Vector3D split_point = centroid_box.centroid();
  std::vector<Primitive *> prims_less;
  std::vector<Primitive *> prims_greater;
  if (split_dim == 0) {
    for (Primitive *p : prims) {
      Vector3D c = p->get_bbox().centroid();
      if (c.x < split_point.x) {
        prims_less.push_back(p);
      } else {
        prims_greater.push_back(p);
      }
    }
  } else if (split_dim == 1) {
    for (Primitive *p : prims) {
      Vector3D c = p->get_bbox().centroid();
      if (c.y < split_point.y) {
        prims_less.push_back(p);
      } else {
        prims_greater.push_back(p);
      }
    }
  } else {
    for (Primitive *p : prims) {
      Vector3D c = p->get_bbox().centroid();
      if (c.z < split_point.z) {
        prims_less.push_back(p);
      } else {
        prims_greater.push_back(p);
      }
    }
  }
  if ((prims_less.size() == 0) || (prims_greater.size() == 0)) {
    node->prims = new vector<Primitive *>(prims);
    return node;
  }
  node->l = construct_bvh(prims_less, max_leaf_size);
  node->r = construct_bvh(prims_greater, max_leaf_size);
  return node;

}


bool BVHAccel::intersect(const Ray& ray, BVHNode *node) const {

  // TODO Part 2, task 3:
  // Implement BVH intersection.
  // Currently, we just naively loop over every primitive.

  double t0 = ray.min_t;
  double t1 = ray.max_t;
  if (!node->bb.intersect(ray, t0, t1)) {
    return false;
  }
  if (node->isLeaf()) {
    for (Primitive *p : *(node->prims)) {
      if (p->intersect(ray)) { 
        return true;
      }
    }
    return false;
  }
  return intersect(ray, node->l) || intersect(ray, node->r);

}

bool BVHAccel::intersect(const Ray& ray, Intersection* i, BVHNode *node) const {

  // TODO Part 2, task 3:
  // Implement BVH intersection.
  // Currently, we just naively loop over every primitive.

  bool hit = false;
  double t0 = ray.min_t;
  double t1 = ray.max_t;
  if (!node->bb.intersect(ray, t0, t1)) {
    return false;
  }
  if (node->isLeaf()) {
    for (Primitive *p : *(node->prims)) {
      if (p->intersect(ray, i)) { 
        total_isects++;
        hit = true;
      }
    }
    return hit;
  }
  bool hit1 = intersect(ray, i, node->l);
  bool hit2 = intersect(ray, i, node->r);
  hit = hit1 || hit2;
  return hit;

}

}  // namespace StaticScene
}  // namespace CGL
