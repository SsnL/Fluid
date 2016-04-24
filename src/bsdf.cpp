#include "bsdf.h"

#include <iostream>
#include <algorithm>
#include <utility>

using std::min;
using std::max;
using std::swap;

namespace CGL {

void make_coord_space(Matrix3x3& o2w, const Vector3D& n) {

    Vector3D z = Vector3D(n.x, n.y, n.z);
    Vector3D h = z;
    if (fabs(h.x) <= fabs(h.y) && fabs(h.x) <= fabs(h.z)) h.x = 1.0;
    else if (fabs(h.y) <= fabs(h.x) && fabs(h.y) <= fabs(h.z)) h.y = 1.0;
    else h.z = 1.0;

    z.normalize();
    Vector3D y = cross(h, z);
    y.normalize();
    Vector3D x = cross(z, y);
    x.normalize();

    o2w[0] = x;
    o2w[1] = y;
    o2w[2] = z;
}

// Diffuse BSDF //

Spectrum DiffuseBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return albedo * (1.0 / PI);
}

Spectrum DiffuseBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *wi = sampler.get_sample(pdf);
  return albedo * (1.0 / PI);
}

// Mirror BSDF //

Spectrum MirrorBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum MirrorBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO Part 5:
  // Implement MirrorBSDF
  reflect(wo, wi);
  *pdf = 1.0;
  return reflectance * (1 / fabs(wi->z));
}

// Glossy BSDF //

/*
Spectrum GlossyBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum GlossyBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *pdf = 1.0f;
  return reflect(wo, wi, reflectance);
}
*/

// Refraction BSDF //

Spectrum RefractionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum RefractionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO Part 5:
  // Implement RefractionBSDF
  return Spectrum();
}

// Glass BSDF //

Spectrum GlassBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum GlassBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO Part 5:
  // Compute Fresnel coefficient and either reflect or refract based on it.
  if (!refract(wo, wi, ior)) {
    reflect(wo, wi);
    *pdf = 1.0;
    return reflectance * (1.0 / fabs(wi->z));
  }
  double R = (ior - 1.0) * (ior - 1.0) / (ior + 1.0) / (ior + 1.0);
  double c = 1.0 - fabs(wi->z);
  R = R + (1.0 - R) * c * c * c * c * c;
  if (coin_flip(R)) {
    reflect(wo, wi);
    *pdf = R;
    return reflectance * (R / fabs(wi->z));
  }
  *pdf = 1.0 - R;
  double nr;
  if (wi->z > 0) {
    nr = ior;
  } else {
    nr = 1.0 / ior;
  }
  return transmittance * ((1.0 - R) * nr * nr / fabs(wi->z));
}

void BSDF::reflect(const Vector3D& wo, Vector3D* wi) {

  // TODO Part 5:
  // Implement reflection of wo about normal (0,0,1) and store result in wi.
  *wi = Vector3D(-wo.x, -wo.y, wo.z);

}

bool BSDF::refract(const Vector3D& wo, Vector3D* wi, float ior) {

  // TODO Part 5:
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
  // ray entering the surface through vacuum.
  double ni, no;
  if (wo.z > 0) {
    no = 1.0;
    ni = ior;
  } else {
    no = ior;
    ni = 1.0;
  }
  double sint2 = (1.0 - wo.z * wo.z) * no * no / ni / ni;
  if (sint2 >= 1.0) {
    return false;
  }
  double cost = sqrt(1.0 - sint2);
  if (wo.z > 0) {
    cost = -cost;
  }
  *wi = Vector3D(-wo.x * no / ni, -wo.y * no / ni, cost);
  wi->normalize();
  return true;

}

// Emission BSDF //

Spectrum EmissionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum EmissionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *pdf = 1.0 / PI;
  *wi  = sampler.get_sample(pdf);
  return Spectrum();
}

} // namespace CGL
