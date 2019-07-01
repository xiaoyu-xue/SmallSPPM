#include "sampling.h"

NAMESPACE_BEGIN

Vec2 ConcentricSampleDisk(const Vec2 &u) {
	// Map uniform Random numbers to $[-1,1]^2$
	Vec2 uOffset = 2.0 * u - Vec2(1, 1);

	// Handle degeneracy at the origin
	if (uOffset.x == 0 && uOffset.y == 0) return Vec2(0, 0);

	// Apply concentric mapping to point
	real theta, r;
	if (std::abs(uOffset.x) > std::abs(uOffset.y)) {
		r = uOffset.x;
		theta = PiOver4 * (uOffset.y / uOffset.x);
	}
	else {
		r = uOffset.y;
		theta = PiOver2 - PiOver4 * (uOffset.x / uOffset.y);
	}
	return r * Vec2(std::cos(theta), std::sin(theta));
}

Vec3 UniformSampleSphere(const Vec2 &u) {
	real z = 1 - 2 * u.x;
	real r = std::sqrt(std::max((real)0, (real)1 - z * z));
	real phi = 2 * PI * u.y;
	return Vec3(r * std::cos(phi), r * std::sin(phi), z);
}

Vec3 CosineSampleHemisphere(const Vec2 &u) {
	Vec2 d = ConcentricSampleDisk(u);
	real z = std::sqrt(std::max((real)0, 1 - d.x * d.x - d.y * d.y));
	return Vec3(d.x, d.y, z);
}

real BalanceHeuristic(int nf, real fPdf, int ng, real gPdf) {
	return (nf * fPdf) / (nf * fPdf + ng * gPdf);
}


real CosineHemispherePdf(real cosTheta) { 
	return cosTheta * INV_PI; 
}

NAMESPACE_END