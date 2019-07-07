#pragma once

#include "linagl.h"
#include "utils.h"
#include "intersection.h"

NAMESPACE_BEGIN

class BSDF {
public:
	BSDF(const Intersection &isect) : n(isect.n), nl(isect.nl) {}
	virtual real Pdf(const Vec3 &wo, const Vec3 &wi) const = 0;
	virtual Vec3 Sample_f(const Vec3 &wo, Vec3 *wi, real *pdf, const Vec3 &rand = Vec3(0, 0, 0)) const = 0;
	virtual Vec3 f(const Vec3 &wo, const Vec3 &wi) const = 0;
	virtual bool IsDelta() const { return false; }
protected:
	const Vec3 n, nl;
};

class DiffuseBSDF : public BSDF {
public:
	DiffuseBSDF(const Intersection &isect, Vec3 r) : BSDF(isect), R(r) {}

	real Pdf(const Vec3 &wo, const Vec3 &wi) const override {
		return std::abs(wi.Dot(nl)) * INV_PI;
	}

	Vec3 Sample_f(const Vec3 &wo, Vec3 *wi, real *pdf, const Vec3 &rand) const override {
		real r1 = 2.f * PI * rand[0], r2 = rand[1];
		real r2s = sqrt(r2);
		Vec3 w = nl, u = ((fabs(w.x) > .1f ? Vec3(0, 1, 0) : Vec3(1, 0, 0)).Cross(w)).Norm();
		Vec3 v = w.Cross(u);
		*wi = (u* cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).Norm();
		*pdf = Pdf(wo, *wi);
		return f(wo, *wi);
	}

	Vec3 f(const Vec3 &wo, const Vec3 &wi) const override {
		return R * INV_PI;
	}
private:
	Vec3 R;
};

class SpecularBSDF : public BSDF {
public:
	SpecularBSDF(const Intersection &isect, Vec3 r = Vec3(1.0, 1.0, 1.0)) : BSDF(isect), R(r) {}

	real Pdf(const Vec3 &wo, const Vec3 &wi) const override {
		return 0.0;
	}

	Vec3 Sample_f(const Vec3 &wo, Vec3 *wi, real *pdf, const Vec3 &rand) const override {
		*wi = (nl * 2.0 * nl.Dot(wo) - wo).Norm();
		*pdf = 1.0;
		real cosTheta = std::abs((*wi).Dot(n));
		return R / cosTheta;
	}

	Vec3 f(const Vec3 &wo, const Vec3 &wi) const override {
		return Vec3();
	}

	bool IsDelta() const override { return true; }
private:
	Vec3 R;
};

class TransmissionBSDF : public BSDF {
public:
	TransmissionBSDF(const Intersection &isect, Vec3 fa = Vec3(1.0, 1.0, 1.0), TransportMode mode = TransportMode::Radiance,
		real eta1 = 1.0, real eta2 = 1.5) :
		BSDF(isect), Fa(fa), nc(eta1), nt(eta2), transportMode(mode) {
	}

	real Pdf(const Vec3 &wo, const Vec3 &wi) const override {
		return 0.0;
	}

	Vec3 Sample_f(const Vec3 &wo, Vec3 *wi, real *pdf, const Vec3 &rand) const override {
		/*
		bool into = (n.Dot(nl) > 0.0);
		real nnt = into ? nc / nt : nt / nc, ddn = (-1 * wo).Dot(nl), cos2t;
		// total internal reflection
		if ((cos2t = 1 - nnt * nnt*(1 - ddn * ddn)) < 0) {
			*wi = (nl * 2.0 * nl.Dot(wo) - wo).Norm();
			real cosTheta = std::abs((*wi).Dot(n));
			*pdf = 1.0;
			return Fa / cosTheta;
			//std::cout << "total internal reflection" << std::endl;
			//return Vec3();
		}
		Vec3 td = ((-1 * wo) * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).Norm();
		real Re = Fresnell(wo, td, n, nl);*/


		real Re = FrDielectric(wo.Dot(n), nc, nt);

		bool entering = wo.Dot(n) > 0;
		real etaI = entering ? nc : nt;
		real etaT = entering ? nt : nc;
		Vec3 nForward = wo.Dot(n) > 0 ? n : -1 * n;
		real P = Re * 0.5f + 0.25f;
		//real P = Re;
		if (rand.z < P) {
			//{
			//	if (debugPixel == 1) {
			//		std::cout << "Reflection" << std::endl;
			//	}
			//}
			/*
			*wi = (nl * 2.0 * nl.Dot(wo) - wo).Norm();
			*pdf = P;
			real cosTheta = std::abs((*wi).Dot(n));
			return Fa * Re / cosTheta;
			*/
			*wi = (2 * wo.Dot(nForward) * nForward - wo).Norm();
			*pdf = P;
			real cosTheta = std::abs((*wi).Dot(nForward));
			return Fa * Re / cosTheta;
		}
		else {
			/*
			*wi = td;
			*pdf = 1 - P;
			real cosTheta = std::abs((*wi).Dot(n));
			return nnt * nnt * Fa * (1.0 - Re) / cosTheta;
			*/
			//{
			//	if (debugPixel == 1) {
			//		std::cout << "Refraction" << std::endl;
			//	}
			//}
			real cosThetaI = wo.Dot(nForward);
			real sinThetaI = std::sqrt(std::max(1 - cosThetaI * cosThetaI, (real)0));
			real sinThetaT = sinThetaI * etaI / etaT;
			if (sinThetaT >= 1.0) {
				*pdf = 1.0;
				*wi = (2 * wo.Dot(nForward) * nForward - wo).Norm();
				return Vec3();
			}
			real cosThetaT = std::sqrt(std::max(1 - sinThetaT * sinThetaT, (real)0));
			*wi = (cosThetaI * nForward - wo).Norm() * sinThetaT - nForward * cosThetaT;
			*pdf = 1 - P;
			real cosTheta = std::abs((*wi).Dot(nForward));
			real eta = etaI / etaT;
			if (transportMode == TransportMode::Radiance) {
				return eta * eta * Fa * (1.f - Re) / cosTheta;
			}
			else {
				return Fa * (1.f - Re) / cosTheta;
			}
		}
	}

	Vec3 f(const Vec3 &wo, const Vec3 &wi) const override {
		return Vec3();
	}

	real Fresnell(const Vec3 &wo, const Vec3 &td, const Vec3 &n, const Vec3 &nl) const {
		bool into = (n.Dot(nl) > 0.0);
		real nnt = into ? nc / nt : nt / nc, ddn = (-1 * wo).Dot(nl);//, cos2t;
		//cos2t = 1 - nnt * nnt*(1 - ddn * ddn);
		real a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : td.Dot(n));
		real Re = R0 + (1 - R0) * c * c* c * c * c;
		return Re;
	}

	real FrDielectric(real cosThetaI, real etaI, real etaT) const {
		cosThetaI = Clamp(cosThetaI, -1.f, 1.f);
		// Potentially swap indices of refraction
		bool entering = cosThetaI > 0.f;
		if (!entering) {
			std::swap(etaI, etaT);
			cosThetaI = std::abs(cosThetaI);
		}

		// Compute _cosThetaT_ using Snell's law
		real sinThetaI = std::sqrt(std::max((real)0, 1 - cosThetaI * cosThetaI));
		real sinThetaT = etaI / etaT * sinThetaI;

		// Handle total internal reflection
		if (sinThetaT >= 1) return 1;
		real cosThetaT = std::sqrt(std::max((real)0, 1 - sinThetaT * sinThetaT));
		real Rparl = ((etaT * cosThetaI) - (etaI * cosThetaT)) /
			((etaT * cosThetaI) + (etaI * cosThetaT));
		real Rperp = ((etaI * cosThetaI) - (etaT * cosThetaT)) /
			((etaI * cosThetaI) + (etaT * cosThetaT));
		return (Rparl * Rparl + Rperp * Rperp) / 2;
	}

	bool IsDelta() const override { return true; }

	bool Refract(const Vec3 &wi, const Vec3 &n, real eta, Vec3 *wt) {
		// Compute $\cos \theta_\roman{t}$ using Snell's law
		real cosThetaI = n.Dot(wi);
		real sin2ThetaI = std::max(real(0), real(1 - cosThetaI * cosThetaI));
		real sin2ThetaT = eta * eta * sin2ThetaI;

		// Handle total internal reflection for transmission
		if (sin2ThetaT >= 1) return false;
		real cosThetaT = std::sqrt(1 - sin2ThetaT);
		*wt = eta * (-1 * wi) + (eta * cosThetaI - cosThetaT) * Vec3(n);
		return true;
	}
private:
	real nc, nt;
	Vec3 Fa;
	TransportMode transportMode;
};


NAMESPACE_END