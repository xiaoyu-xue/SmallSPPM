#include "Homogeneous.h"
#include "core/Intersection.h"
#include "core/Scene.h"
#include <algorithm>
GYT_NAMESPACE_BEGIN

Vec3 HomogeneousMedium::Tr(const Ray& ray, StateSequence& rand) const {
    return Exp(-sigma_t * std::min(ray.m_tMax * ray.mDir.Length(), MaxReal));
}


Vec3 HomogeneousMedium::Sample(const Ray& ray, StateSequence& rand, MemoryPool& arena, MediumIntersection* mi) const {
    // Sample a channel and distance along the ray
    int channel = std::min((int)(rand() * 3), 2);
    //int channel = 0;
    real dist = -std::log(1 - rand()) / sigma_t[channel];
    //real dist = -std::log(1 - std::max(0.f, rand() - 0.1f)) / sigma_t[channel];
    real t = std::min(dist / ray.mDir.Length(), ray.m_tMax);
    bool sampledMedium = t < ray.m_tMax;
    if (sampledMedium)
        *mi = MediumIntersection(ray(t), -ray.mDir, this, MEMORY_POOL_ALLOC(arena, HenyeyGreenstein)(g));

    // Compute the transmittance and sampling density
    Vec3 Tr = Exp(-sigma_t * std::min(t, MaxReal) * ray.mDir.Length());

    // Return weighting factor for scattering from homogeneous medium
    Vec3 density = sampledMedium ? (sigma_t * Tr) : Tr;
    real pdf = 0;
    for (int i = 0; i < 3; ++i) pdf += density[i];
    pdf *= 1 / (real)(3);
    //pdf = density[0];
    if (pdf == 0) {
        pdf = 1;
    }
    return sampledMedium ? (Tr * sigma_s / pdf) : (Tr / pdf);
}

Vec3 HomogeneousMedium::EquiAngularSampling(
    const Ray& ray, StateSequence& rand, MemoryPool& arena,
    const Intersection& lightPoint, MediumIntersection* mi) const {
    //real rrProb = Exp(-sigma_t * std::min(ray.tMax * ray.d.Length(), MaxReal)).x;
    //real rrProb = Tr(ray, rand).x;
    //if (rand() > rrProb) return 1.f / (1 - rrProb);
    real tMax = 3;
    //real tMax = std::min((real)10, ray.tMax);


    Vec3 lightPos = lightPoint.mPos;
    real delta = Dot(lightPos - ray.mOrig, ray.mDir);
    real D = (ray.mOrig + delta * ray.mDir - lightPos).Length();
    real theta_a = std::atan2(0.f - delta, D);
    real theta_b = std::atan2(tMax - delta, D);
    real u = rand();
    real tau = D * std::tan((1 - u) * theta_a + u * theta_b);
    real t = delta + tau;
    bool sampleMedia = t < ray.m_tMax;
    Vec3 Tr = sampleMedia ? Exp(-sigma_t * t) : this->Tr(ray, rand);
    real pdf = sampleMedia ?
        D / ((theta_b - theta_a) * (D * D + tau * tau))
        : 1 - (std::atan2(std::min(tMax, ray.m_tMax) - delta, D) - theta_a) / (theta_b - theta_a);
    if(sampleMedia) *mi = MediumIntersection(ray(t), -ray.mDir, this, MEMORY_POOL_ALLOC(arena, HenyeyGreenstein)(g));

    return sampleMedia ? Tr * sigma_s / pdf : this->Tr(ray, rand) / pdf;
}
GYT_NAMESPACE_END