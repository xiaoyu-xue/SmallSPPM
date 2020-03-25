#include "homogeneous.h"
#include <algorithm>
#include "intersection.h"
#include "scene.h"

NAMESPACE_BEGIN

Vec3 HomogeneousMedium::Tr(const Ray& ray, StateSequence& rand) const {
    return Exp(-sigma_t * std::min(ray.tMax * ray.d.Length(), MaxReal));
}


Vec3 HomogeneousMedium::Sample(const Ray& ray, StateSequence& rand, MemoryArena& arena, MediumIntersection* mi) const {
    // Sample a channel and distance along the ray
    int channel = std::min((int)(rand() * 3), 2);
    //int channel = 0;
    real dist = -std::log(1 - rand()) / sigma_t[channel];
    //real dist = -std::log(1 - std::max(0.f, rand() - 0.1f)) / sigma_t[channel];
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
    //pdf = density[0];
    if (pdf == 0) {
        pdf = 1;
    }
    return sampledMedium ? (Tr * sigma_s / pdf) : (Tr / pdf);
}

Vec3 HomogeneousMedium::EquiAngularSampling(
    const Ray& ray, StateSequence& rand, MemoryArena& arena,
    const Intersection& lightPoint, MediumIntersection* mi) const {
    //real rrProb = Exp(-sigma_t * std::min(ray.tMax * ray.d.Length(), MaxReal)).x;
    ////real rrProb = 0.8;
    //if (rand() > rrProb) return 1.f / (1 - rrProb);
    real tMax = std::min((real)10, ray.tMax);

    //Vec3 c = lightPoint.hit, a = ray.o, b = ray(tMax);
    //Vec3 ac = c - a, ca = a - c, cb = b - c, ab = b - a;
    //Vec3 proj = Dot(ac, ab) * ab / ab.Length2();
    //Vec3 h = proj - ac;
    //real D = h.Length();
    //real cosTheta_b = Dot(h.Norm(), cb.Norm());
    //real theta_b = std::acos(cosTheta_b);
    //real theta_a = Dot(h.Norm(), ca.Norm());
    //theta_b *= -Sgn(Dot(cb, ray.d));
    //theta_a *= -Sgn(Dot(ca, ray.d));
    //real u = rand();
    //real tau = D * std::tan((1 - u) * theta_a + u * theta_b);
    //real t = std::sqrt(std::max(real(0), ca.Length2() - D * D)) + tau;
    //real pdf = D / ((theta_b - theta_a) * (D * D + t * t));
    Vec3 lightPos = lightPoint.hit;
    real delta = Dot(lightPos - ray.o, ray.d);
    real D = (ray.o + delta * ray.d - lightPos).Length();
    real theta_a = std::atan2(0.f - delta, D);
    real theta_b = std::atan2(tMax - delta, D);
    real u = rand();
    real tau = D * std::tan((1 - u) * theta_a + u * theta_b);
    real pdf = D / ((theta_b - theta_a) * (D * D + tau * tau));
    real t = delta + tau;
    Vec3 Tr = Exp(-sigma_t * t);
    *mi = MediumIntersection(ray(t), -ray.d, this, ARENA_ALLOC(arena, HenyeyGreenstein)(g));

    //std::cout << t << " " << tMax << std::endl;
    if (pdf < 0) std::cout << "pdf " << pdf << std::endl;
    return Tr * sigma_s / pdf;
}
NAMESPACE_END