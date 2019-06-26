#include "sampling.h"

NAMESPACE_BEGIN

Vec ConcentricSampleDisk(const Vec &u) {
	// Map uniform Random numbers to $[-1,1]^2$
	Vec uOffset = 2.f * u - Vec(1, 1);

	// Handle degeneracy at the origin
	if (uOffset.x == 0 && uOffset.y == 0) return Vec(0, 0);

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
	return r * Vec(std::cos(theta), std::sin(theta));
}

Vec UniformSampleSphere(const Vec &u) {
	real z = 1 - 2 * u.x;
	real r = std::sqrt(std::max((real)0, (real)1 - z * z));
	real phi = 2 * PI * u.y;
	return Vec(r * std::cos(phi), r * std::sin(phi), z);
}

Vec CosineSampleHemisphere(const Vec &u) {
	Vec d = ConcentricSampleDisk(u);
	real z = std::sqrt(std::max((real)0, 1 - d.x * d.x - d.y * d.y));
	return Vec(d.x, d.y, z);
}

real BalanceHeuristic(int nf, real fPdf, int ng, real gPdf) {
	return (nf * fPdf) / (nf * fPdf + ng * gPdf);
}


real CosineHemispherePdf(real cosTheta) { 
	return cosTheta * INV_PI; 
}

NAMESPACE_END