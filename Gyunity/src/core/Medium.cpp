#include "Medium.h"
#include "math/GeometryUtils.h"
#include <algorithm>

GYT_NAMESPACE_BEGIN

// HenyeyGreenstein Method Definitions
real HenyeyGreenstein::Sample_p(const Vec3& wo, Vec3* wi, const Vec2& u) const {
    // Compute $\cos \theta$ for Henyey--Greenstein sample
    real cosTheta;
    if (std::abs(g) < 1e-3)
        cosTheta = 1 - 2 * u[0];
    else {
        real sqrTerm = (1 - g * g) / (1 - g + 2 * g * u[0]);
        cosTheta = (1 + g * g - sqrTerm * sqrTerm) / (2 * g);
    }

    // Compute direction _wi_ for Henyey--Greenstein sample
    real sinTheta = std::sqrt(std::max((real)0, 1 - cosTheta * cosTheta));
    real phi = 2 * PI * u[1];
    Vec3 v1, v2;
    CoordinateSystem(wo, &v1, &v2);
    *wi = SphericalDirection(sinTheta, cosTheta, phi, v1, v2, -wo);
    return PhaseHG(-cosTheta, g);
}

real HenyeyGreenstein::p(const Vec3& wo, const Vec3& wi) const {
    return PhaseHG(Dot(wo, wi), g);
}

GYT_NAMESPACE_END