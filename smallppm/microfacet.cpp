#include "microfacet.h"
#include "bsdf.h"

NAMESPACE_BEGIN

static void TrowbridgeReitzSample11(real cosTheta, real U1, real U2,
	real* slope_x, real* slope_y) {
	// special case (normal incidence)
	if (cosTheta > .9999) {
		real r = sqrt(U1 / (1 - U1));
		real phi = 6.28318530718 * U2;
		*slope_x = r * cos(phi);
		*slope_y = r * sin(phi);
		return;
	}

	real sinTheta =
		std::sqrt(std::max((real)0, (real)1 - cosTheta * cosTheta));
	real tanTheta = sinTheta / cosTheta;
	real a = 1 / tanTheta;
	real G1 = 2 / (1 + std::sqrt(1.f + 1.f / (a * a)));

	// sample slope_x
	real A = 2 * U1 / G1 - 1;
	real tmp = 1.f / (A * A - 1.f);
	if (tmp > 1e10) tmp = 1e10;
	real B = tanTheta;
	real D = std::sqrt(
		std::max(real(B * B * tmp * tmp - (A * A - B * B) * tmp), real(0)));
	real slope_x_1 = B * tmp - D;
	real slope_x_2 = B * tmp + D;
	*slope_x = (A < 0 || slope_x_2 > 1.f / tanTheta) ? slope_x_1 : slope_x_2;

	// sample slope_y
	real S;
	if (U2 > 0.5f) {
		S = 1.f;
		U2 = 2.f * (U2 - .5f);
	}
	else {
		S = -1.f;
		U2 = 2.f * (.5f - U2);
	}
	real z =
		(U2 * (U2 * (U2 * 0.27385f - 0.73369f) + 0.46341f)) /
		(U2 * (U2 * (U2 * 0.093073f + 0.309420f) - 1.000000f) + 0.597999f);
	*slope_y = S * z * std::sqrt(1.f + *slope_x * *slope_x);

}

static Vec3 TrowbridgeReitzSample(const Vec3& wi, real alpha_x,
	real alpha_y, real U1, real U2) {
	// 1. stretch wi
	Vec3 wiStretched =
		Vec3(alpha_x * wi.x, alpha_y * wi.y, wi.z).Norm();

	// 2. simulate P22_{wi}(x_slope, y_slope, 1, 1)
	real slope_x, slope_y;
	TrowbridgeReitzSample11(BSDFCoordinate::CosTheta(wiStretched), U1, U2, &slope_x, &slope_y);

	// 3. rotate
	real tmp = BSDFCoordinate::CosPhi(wiStretched) * slope_x - BSDFCoordinate::SinPhi(wiStretched) * slope_y;
	slope_y = BSDFCoordinate::SinPhi(wiStretched) * slope_x + BSDFCoordinate::CosPhi(wiStretched) * slope_y;
	slope_x = tmp;

	// 4. unstretch
	slope_x = alpha_x * slope_x;
	slope_y = alpha_y * slope_y;

	// 5. compute normal
	return Vec3(-slope_x, -slope_y, 1.).Norm();
}


real MicrofacetDistribution::Pdf(const Vec3& wo, const Vec3& wh) const {
    if (sampleVisibleArea)
        return D(wh) * G1(wo, wh) * std::abs(Dot(wo, wh)) / BSDFCoordinate::AbsCosTheta(wo);
    else
        return D(wh) * BSDFCoordinate::AbsCosTheta(wh);
}

real GGXDistribution::Lambda(const Vec3& w) const {
	real tanTheta = BSDFCoordinate::TanTheta(w);
	if (std::isinf(std::abs(tanTheta))) return 0.f;
	real inv_a2;
	if (alphax == alphay) {
		inv_a2 = alphax * alphax * tanTheta * tanTheta;
	}
	else {
		real cos2Phi = BSDFCoordinate::Cos2Phi(w);
		real sin2Phi = BSDFCoordinate::Sin2Phi(w);
		real alpha_o = std::sqrt(cos2Phi * alphax * alphax + sin2Phi * alphay * alphay);
		inv_a2 = alpha_o * alpha_o * tanTheta * tanTheta;
	}
	return (-1 + std::sqrt(1 + inv_a2)) / 2;
}

real GGXDistribution::D(const Vec3& wh) const {
	if (alphax == alphay) {
		real alpha = alphax;
		const real alpha2 = alpha * alpha;
		const real cos2ThetaM = BSDFCoordinate::Cos2Theta(wh);
		const real cos4ThetaM = cos2ThetaM * cos2ThetaM;
		const real tan2ThetaM = BSDFCoordinate::Tan2Theta(wh);

		if (std::isinf(tan2ThetaM)) return 0.;
		const real root = alpha2 + tan2ThetaM;

		return alpha2 / (PI * cos4ThetaM * root * root);
	}
	else {
		const real cos2PhiM = BSDFCoordinate::Cos2Phi(wh);
		const real sin2PhiM = BSDFCoordinate::Sin2Phi(wh);
		const real tan2ThetaM = BSDFCoordinate::Tan2Theta(wh);
		const real cosThetaM = BSDFCoordinate::CosTheta(wh);
		const real cos2ThetaM = BSDFCoordinate::Cos2Theta(wh);

		if (std::isinf(tan2ThetaM)) return 0.;
		const real alphax2 = alphax * alphax;
		const real alphay2 = alphay * alphay;
		const real root = 1 + tan2ThetaM * (cos2PhiM / alphax2 + sin2PhiM / alphay2);
		return 1 / (PI * alphax * alphay * cos2ThetaM * cos2ThetaM * root * root);
	}

}

real GGXDistribution::G1(const Vec3& v, const Vec3& wh) const {

	if (alphax == alphay) {
		real alpha = alphax;
		const real tanTheta = std::abs(BSDFCoordinate::TanTheta(v));
		if (tanTheta == 0.0f)
			return 1.0f;

		if (Dot(v, wh) * BSDFCoordinate::CosTheta(v) <= 0)
			return 0.0f;

		const real root = alpha * tanTheta;
		return 2.0f / (1.0f + std::sqrt(1.0f + root * root));
	}
	else {
		return 1 / (1 + Lambda(v));
	}

}

real GGXDistribution::G(const Vec3& wo, const Vec3& wi, const Vec3& wh) const {
    //return G1(wo, wh) * G1(wi, wh);
	return 1 / (1 + Lambda(wo) + Lambda(wi));
}

Vec3 GGXDistribution::Sample_wh(const Vec3& wo, const Vec2& u) const {
	if (!sampleVisibleArea) {
		if (alphax == alphay) {
			real alpha = alphax;
			real alpha2 = alpha * alpha;
			real tan2Theta = alpha2 * u[0] / (1 - u[0]);
			real cosTheta = 1.f / std::sqrt(1 + tan2Theta);
			real cos2Theta = cosTheta * cosTheta;
			real sinTheta = std::sqrt(std::max(real(0), 1 - cos2Theta));
			real phi = 2 * PI * u[1];
			Vec3 dir = SphericalDirection(sinTheta, cosTheta, phi);
			return dir;
		}
		else {
			real tanPhi = alphay / alphax * std::tan(2 * PI * u[1]);
			real phi;
			if (u[1] >= 0 && u[1] <= 0.25) {
				phi = std::atan(tanPhi);
			}
			else if (u[1] > 0.25 && u[1] < 0.75) {
				phi = std::atan(tanPhi) + PI;
			}
			else {
				phi = std::atan(tanPhi) + 2 * PI;
			}
			real cosPhi = std::cos(phi);
			real cos2Phi = cosPhi * cosPhi;
			real sinPhi = std::sin(phi);
			real sin2Phi = sinPhi * sinPhi;
			real alphax2 = alphax * alphax;
			real alphay2 = alphay * alphay;
			real tanTheta = u[0] / ((1 - u[0]) * (cos2Phi / alphax2 + sin2Phi / alphay2));
			real theta = std::atan(tanTheta);
			real cosTheta = std::cos(theta);
			real sinTheta = std::cos(theta);
			Vec3 dir = SphericalDirection(sinTheta, cosTheta, phi);
			return dir;
		}
	}
	else {
		//return SampleGGXVNDF(wo, u);
		bool flip = wo.z < 0;
		Vec3 wh = TrowbridgeReitzSample(flip ? -wo : wo, alphax, alphay, u[0], u[1]);
		if (flip) wh = -wh;
		return wh;
	}
	

}

//https://hal.archives-ouvertes.fr/hal-01509746
Vec3 GGXDistribution::SampleGGXVNDF(const Vec3& v_, const Vec2& u) const {
	//stretch view
	Vec3 v = Vec3(alphax * v_.x, alphay * v_.y, v_.z).Norm();

	//orthonormal basis
	Vec3 T1 = (v.z < 0.9999) ? Cross(v, Vec3(0, 0, 1)).Norm() : Vec3(1, 0, 0);
	Vec3 T2 = Cross(T1, v);

	//sample point with polar coordinates
	real a = 1 / (1.f + v.z);
	real r = std::sqrt(u[0]);
	real phi = (u[1] < a) ? u[1] / a * PI : PI + (u[1] - a) / (1.f - a) * PI;
	real P1 = r * std::cos(phi);
	real P2 = r * std::sin(phi) * ((u[1] < a) ? 1.0f : v.z);

	//compute normal
	Vec3 dir = P1 * T1 + P2 * T2 + std::sqrt(std::max(real(0), 1 - P1 * P1 - P2 * P2)) * v;

	//unstretch
	dir = Vec3(alphax * dir.x, alphay * dir.y, std::max(real(0), dir.z)).Norm();
	return dir;
}

NAMESPACE_END