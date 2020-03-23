#include "homogeneous.h"
#include <algorithm>
#include "intersection.h"

NAMESPACE_BEGIN

Vec3 HomogeneousMedium::Tr(const Ray& ray, StateSequence& rand) const {
    return Exp(-sigma_t * std::min(ray.tMax * ray.d.Length(), MaxReal));
}


Vec3 HomogeneousMedium::Sample(const Ray& ray, StateSequence& rand, MemoryArena& arena, MediumIntersection* mi) const {
    // Sample a channel and distance along the ray
    int channel = (int)(rand() * 3);
    real dist = -std::log(1 - rand()) / sigma_t[channel];
    real t = std::min(dist / ray.d.Length(), ray.tMax);
    bool sampledMedium = t < ray.tMax;
    if (sampledMedium)
        *mi = MediumIntersection(ray(t), -ray.d, this, ARENA_ALLOC(arena, HenyeyGreenstein)(g));

    // Compute the transmittance and sampling density
    Vec3 Tr = Exp(-sigma_t * std::min(t, MaxReal) * ray.d.Length());

    // Return weighting factor for scattering from homogeneous medium
    Vec3 density = sampledMedium ? (sigma_t * Tr) : Tr;
    real pdf = 0;
    for (int i = 0; i < 3; ++i) pdf += density[i];
    pdf *= 1 / (real)(3);
    if (pdf == 0) {
        pdf = 1;
    }
    return sampledMedium ? (Tr * sigma_s / pdf) : (Tr / pdf);
}
NAMESPACE_END