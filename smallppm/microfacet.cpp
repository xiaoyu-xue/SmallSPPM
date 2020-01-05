#include "microfacet.h"
#include "bsdf.h"

NAMESPACE_BEGIN

real MicrofacetDistribution::Pdf(const Vec3& wo, const Vec3& wh) const {
    if (sampleVisibleArea)
        return D(wh) * G1(wo, wh) * std::abs(Dot(wo, wh)) / BSDFCoordinate::AbsCosTheta(wo);
    else
        return D(wh) * BSDFCoordinate::AbsCosTheta(wh);
}

real GGXDistribution::D(const Vec3& wh) const {
    const real alpha2 = alpha * alpha;
    const real cos2ThetaM = BSDFCoordinate::Cos2Theta(wh);
    const real cos4ThetaM = cos2ThetaM * cos2ThetaM;
    const real tan2ThetaM = BSDFCoordinate::Tan2Theta(wh);

    if (std::isinf(tan2ThetaM)) return 0.;
    const real root = alpha2 + tan2ThetaM;

    return alpha2 / (PI * cos4ThetaM * root * root);
}

real GGXDistribution::G1(const Vec3& v, const Vec3& wh) const {

    const real tanTheta = std::abs(BSDFCoordinate::TanTheta(v));
    if (tanTheta == 0.0f)
        return 1.0f;

    if (Dot(v, wh) * BSDFCoordinate::CosTheta(v) <= 0)
        return 0.0f;

    const real root = alpha * tanTheta;
    return 2.0f / (1.0f + std::sqrt(1.0f + root * root));
}

real GGXDistribution::G(const Vec3& wo, const Vec3& wi, const Vec3& wh) const {
    return G1(wo, wh) * G1(wi, wh);
}

Vec3 GGXDistribution::Sample_wh(const Vec3& wo, const Vec2& u) const {
    real alpha2 = alpha * alpha;
    real tan2Theta = alpha2 * u[0] / (1 - u[0]);
    real cosTheta = 1.f / std::sqrt(1 + tan2Theta);
    real cos2Theta = cosTheta * cosTheta;
    real sinTheta = std::sqrt(std::max(real(0), 1 - cos2Theta));
    real phi = 2 * PI * u[1];
    Vec3 dir = SphericalDirection(sinTheta, cosTheta, phi);
    return dir;
}

NAMESPACE_END