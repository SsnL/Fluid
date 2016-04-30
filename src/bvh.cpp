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
  if (prims) delete prims;
}

BBox BVHAccel::get_bbox() const {
  return root->bb;
}

void BVHAccel::draw(BVHNode *node, const Color& c) const {
  if (node->isLeaf()) {
    for (size_t i = node->l_ind; i < node->r_ind; i++)
      prims->at(i)->draw(c);
  } else {
    draw(node->l, c);
    draw(node->r, c);
  }
}

void BVHAccel::drawOutline(BVHNode *node, const Color& c) const {
  if (node->isLeaf()) {
    for (size_t i = node->l_ind; i < node->r_ind; i++)
      prims->at(i)->drawOutline(c);
  } else {
    drawOutline(node->l, c);
    drawOutline(node->r, c);
  }
}

BVHNode *BVHAccel::construct_bvh_in_range(
  size_t max_leaf_size,
  size_t l, size_t r, BBox &bbox
) {
  size_t n = r - l;

  // If reach threshold size, build a leafs
  if (n <= max_leaf_size) {
    BVHNode *node = new BVHNode(bbox);
    node->l_ind = l;
    node->r_ind = r;
    return node;
  }

  vector<Primitive*>::iterator it_l = prims->begin() + l,
    it_r = prims->begin() + r;

  // Find longest dimension of bounding box
  size_t axis;
  Vector3D extent = bbox.extent;
  if (extent[0] > extent[1] && extent[0] > extent[2]) {
    axis = 0;
  } else if (extent[1] > extent[0] && extent[1] > extent[2]) {
    axis = 1;
  } else {
    axis = 2;
  }

  int i;

  // Sort by comparator on given axis
  sort(it_l, it_r, PrimitiveCentroidComparator(axis));

  // Build bboxs of right child candidates
  // r_bboxes[i] = bbox(prims[l + i:r])
  vector<BBox> r_bboxes(n);
  BBox r_bbox;
  for (i = r - 1; ; i--) {
    r_bbox.expand(prims->at(i)->get_bbox());
    r_bboxes[i - l] = r_bbox;
    if (i == l) {
      break;
    }
  }

  // Find optimum
  // At iteration i, split as [l:i), [i:r)
  // l_bbox = bbox(prims[l:l + i])
  BBox l_bbox = prims->at(l)->get_bbox(), l_bbox_min, r_bbox_min;
  double sa_l, sa_r, heuristic, min = INF_D;
  size_t n_l, min_i;
  for (i = l + 1, n_l = 1; i < r; i++, n_l++) {
    sa_l = l_bbox.surface_area();
    sa_r = r_bboxes[n_l].surface_area();
    heuristic = sa_l * ((double) n_l) + sa_r * ((double) (n - n_l));
    if (heuristic < min) {
      min = heuristic;
      min_i = i;
      l_bbox_min = l_bbox;
      r_bbox_min = r_bboxes[n_l];
    }
    l_bbox.expand(prims->at(i)->get_bbox());
  }

  // Build node
  BVHNode *node = new BVHNode(bbox);
  // size_t min_i = l + (n >> 1);
  // BBox l_bbox_min, r_bbox_min;
  // for (i = l; i < r; i++) {
  //   if (i < min_i) {
  //     l_bbox_min.expand(prims->at(i)->get_bbox());
  //   } else {
  //     r_bbox_min.expand(prims->at(i)->get_bbox());
  //   }
  // }
  node->l = construct_bvh_in_range(max_leaf_size, l, min_i, l_bbox_min);
  node->r = construct_bvh_in_range(max_leaf_size, min_i, r, r_bbox_min);
  return node;
}


BVHNode *BVHAccel::construct_bvh(
  const vector<Primitive*>& prims, size_t max_leaf_size
) {
  BBox bbox;
  for (Primitive *p : prims) {
    bbox.expand(p->get_bbox());
  }
  this->prims = new vector<Primitive *>(prims);
  return construct_bvh_in_range(max_leaf_size, 0, prims.size(), bbox);
}


bool BVHAccel::intersect(const Ray& ray, BVHNode *node) const {
  double t0, t1;
  if (!node->bb.intersect(ray, t0, t1)) {
    return false;
  }
  if (t1 < ray.min_t || t0 > ray.max_t) {
    return false;
  }

  if (node->isLeaf()) {
    for (size_t i = node->l_ind; i < node->r_ind; i++) {
      total_isects++;
      if (prims->at(i)->intersect(ray)) {
        return true;
      }
    }
    return false;
  }

  return intersect(ray, node->l) || intersect(ray, node->r);

}

bool BVHAccel::intersect(const Ray& ray, Intersection* isect, BVHNode *node) const {
  double t0, t1;
  if (!node->bb.intersect(ray, t0, t1)) {
    return false;
  }
  if (t1 < ray.min_t || t0 > ray.max_t) {
    return false;
  }

  bool hit = false;

  if (node->isLeaf()) {
    for (size_t i = node->l_ind; i < node->r_ind; i++) {
      total_isects++;
      if (prims->at(i)->intersect(ray, isect)) {
        hit = true;
      }
    }
    return hit;
  }

  hit = intersect(ray, isect, node->l);
  hit |= intersect(ray, isect, node->r);

  return hit;


}

}  // namespace StaticScene
}  // namespace CGL
