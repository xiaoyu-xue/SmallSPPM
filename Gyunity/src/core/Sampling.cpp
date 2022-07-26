#include "Sampling.h"

GYT_NAMESPACE_BEGIN

Vec2 ConcentricSampleDisk(const Vec2 &u) {
	// Map uniform Random numbers to $[-1,1]^2$
	Vec2 uOffset = 2.0 * u - Vec2(1, 1);

	// Handle degeneracy at the origin
	if (uOffset.x == 0 && uOffset.y == 0) return Vec2(0, 0);

	// Apply concentric mapping to point
	real theta, r;
	if (std::abs(uOffset.x) > std::abs(uOffset.y)) {
		r = uOffset.x;
		theta = PI_Over4 * (uOffset.y / uOffset.x);
	}
	else {
		r = uOffset.y;
		theta = PI_Over2 - PI_Over4 * (uOffset.x / uOffset.y);
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

Vec2 UniformSampleTriangle(const Vec2& u) {
	real su0 = std::sqrt(u[0]);
	return Vec2(1 - su0, u[1] * su0);
}

real BalanceHeuristic(int nf, real fPdf, int ng, real gPdf) {
	return (nf * fPdf) / (nf * fPdf + ng * gPdf);
}


real CosineHemispherePdf(real cosTheta) { 
	return cosTheta * INV_PI; 
}

real UniformSpherePdf() {
	return INV_4PI;
}

Vec3 UniformSampleHemisphere(const Vec2& u)
{
	real x = std::cos(2 * PI * u.y) * std::sqrt(1 - u.x * u.x);
	real y = std::sin(2 * PI * u.y) * std::sqrt(1 - u.x * u.x);
	real z = u.x;
	return Vec3(x, y, z);
}

real UniformSampleHemispherePdf()
{
	return INV_2PI;
}

Vec3 UniformSampleCone(const Vec3& u, real cosThetaMax)
{
	real cosTheta = ((real)1 - u[0]) + u[0] * cosThetaMax;
	real sinTheta = std::sqrt((real)1 - cosTheta * cosTheta);
	real phi = u[1] * 2 * PI;
	return Vec3(std::cos(phi) * sinTheta, std::sin(phi) * sinTheta, cosTheta);
}

Vec3 UniformSampleCone(const Vec3& u, real cosThetaMax, const Vec3& x, const Vec3& y, const Vec3& z)
{
	real cosTheta = Lerp(u[0], cosThetaMax, (real)1.0);
	real sinTheta = std::sqrt((real)1. - cosTheta * cosTheta);
	real phi = u[1] * 2 * PI;
	return std::cos(phi) * sinTheta * x + std::sin(phi) * sinTheta * y + cosTheta * z;
}

real UniformConePdf(real cosThetaMax)
{
	return 1 / (2 * PI * (1 - cosThetaMax));
}

GYT_NAMESPACE_END