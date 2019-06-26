#pragma once

#include "linagl.h"
#include "utils.h"
#include "intersection.h"

NAMESPACE_BEGIN

class BSDF {
public:
	BSDF(const Intersection &isect) : n(isect.n), nl(isect.nl) {}
	virtual real Pdf(const Vec &wo, const Vec &wi) const = 0;
	virtual Vec Sample_f(const Vec &wo, Vec *wi, real *pdf, Vec rand = Vec(0, 0, 0)) const = 0;
	virtual Vec f(const Vec &wo, const Vec &wi) const { return Vec(0, 0, 0); }
	virtual bool IsDelta() const { return false; }
protected:
	const Vec n, nl;
};

class DiffuseBSDF : public BSDF {
public:
	DiffuseBSDF(const Intersection &isect, Vec r) : BSDF(isect), R(r) {}

	real Pdf(const Vec &wo, const Vec &wi) const {
		return std::abs(wi.dot(nl)) * INV_PI;
	}

	Vec Sample_f(const Vec &wo, Vec *wi, real *pdf, Vec rand) const {
		real r1 = 2. * PI * rand[0], r2 = rand[1];
		real r2s = sqrt(r2);
		Vec w = nl, u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).norm();
		Vec v = w % u;
		*wi = (u* cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();
		*pdf = Pdf(wo, *wi);
		return f(wo, *wi);
	}

	Vec f(const Vec &wo, const Vec &wi) const {
		return R * INV_PI;
	}
private:
	Vec R;
};

class SpecularBSDF : public BSDF {
public:
	SpecularBSDF(const Intersection &isect, Vec r = Vec(1.0, 1.0, 1.0)) : BSDF(isect), R(r) {}

	real Pdf(const Vec &wo, const Vec &wi) const {
		return 0.0;
	}

	Vec Sample_f(const Vec &wo, Vec *wi, real *pdf, Vec rand) const {
		*wi = (nl * 2.0 * nl.dot(wo) - wo).norm();
		*pdf = 1.0;
		real cosTheta = std::abs((*wi).dot(n));
		return R / cosTheta;
	}

	Vec f(const Vec &wo, const Vec &wi) const {
		return Vec();
	}

	bool IsDelta() const { return true; }
private:
	Vec R;
};

class TransmissionBSDF : public BSDF {
public:
	TransmissionBSDF(const Intersection &isect, Vec fa = Vec(1.0, 1.0, 1.0), TransportMode mode = TransportMode::Radiance,
		real eta1 = 1.0, real eta2 = 1.5) :
		BSDF(isect), Fa(fa), nc(eta1), nt(eta2), transportMode(mode) {
	}

	real Pdf(const Vec &wo, const Vec &wi) const {
		return 0.0;
	}

	Vec Sample_f(const Vec &wo, Vec *wi, real *pdf, Vec rand) const {
		/*
		bool into = (n.dot(nl) > 0.0);
		real nnt = into ? nc / nt : nt / nc, ddn = (-1 * wo).dot(nl), cos2t;
		// total internal reflection
		if ((cos2t = 1 - nnt * nnt*(1 - ddn * ddn)) < 0) {
			*wi = (nl * 2.0 * nl.dot(wo) - wo).norm();
			real cosTheta = std::abs((*wi).dot(n));
			*pdf = 1.0;
			return Fa / cosTheta;
			//std::cout << "total internal reflection" << std::endl;
			//return Vec();
		}
		Vec td = ((-1 * wo) * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
		real Re = Fresnell(wo, td, n, nl);*/


		real Re = FrDielectric(wo.dot(n), nc, nt);

		bool entering = wo.dot(n) > 0;
		real etaI = entering ? nc : nt;
		real etaT = entering ? nt : nc;
		Vec nForward = wo.dot(n) > 0 ? n : -1 * n;
		real P = Re * 0.5 + 0.25;
		if (rand.z < P) {
			/*
			*wi = (nl * 2.0 * nl.dot(wo) - wo).norm();
			*pdf = P;
			real cosTheta = std::abs((*wi).dot(n));
			return Fa * Re / cosTheta;
			*/
			*wi = (2 * wo.dot(nForward) * nForward - wo).norm();
			*pdf = P;
			real cosTheta = std::abs((*wi).dot(nForward));
			return Fa * Re / cosTheta;
		}
		else {
			/*
			*wi = td;
			*pdf = 1 - P;
			real cosTheta = std::abs((*wi).dot(n));
			return nnt * nnt * Fa * (1.0 - Re) / cosTheta;
			*/
			real cosThetaI = wo.dot(nForward);
			real sinThetaI = std::sqrt(std::max(1 - cosThetaI * cosThetaI, (real)0));
			real sinThetaT = sinThetaI * etaI / etaT;
			if (sinThetaT >= 1.0) {
				*pdf = 1.0;
				*wi = (2 * wo.dot(nForward) * nForward - wo).norm();
				return Vec();
			}
			real cosThetaT = std::sqrt(std::max(1 - sinThetaT * sinThetaT, (real)0));
			*wi = (cosThetaI * nForward - wo).norm() * sinThetaT - nForward * cosThetaT;
			*pdf = 1 - P;
			real cosTheta = std::abs((*wi).dot(nForward));
			real eta = etaI / etaT;
			if (transportMode == TransportMode::Radiance) {
				return eta * eta * Fa * (1.0 - Re) / cosTheta;
			}
			else {
				return Fa * (1.0 - Re) / cosTheta;
			}
		}
	}

	Vec f(const Vec &wo, const Vec &wi) const {
		return Vec();
	}

	real Fresnell(const Vec &wo, const Vec &td, const Vec &n, const Vec &nl) const {
		bool into = (n.dot(nl) > 0.0);
		real nnt = into ? nc / nt : nt / nc, ddn = (-1 * wo).dot(nl);//, cos2t;
		//cos2t = 1 - nnt * nnt*(1 - ddn * ddn);
		real a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : td.dot(n));
		real Re = R0 + (1 - R0) * c * c* c * c * c;
		return Re;
	}

	real FrDielectric(real cosThetaI, real etaI, real etaT) const {
		cosThetaI = Clamp(cosThetaI, -1, 1);
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

	bool IsDelta() const { return true; }

	bool Refract(const Vec &wi, const Vec &n, real eta, Vec *wt) {
		// Compute $\cos \theta_\roman{t}$ using Snell's law
		real cosThetaI = n.dot(wi);
		real sin2ThetaI = std::max(real(0), real(1 - cosThetaI * cosThetaI));
		real sin2ThetaT = eta * eta * sin2ThetaI;

		// Handle total internal reflection for transmission
		if (sin2ThetaT >= 1) return false;
		real cosThetaT = std::sqrt(1 - sin2ThetaT);
		*wt = eta * (-1 * wi) + (eta * cosThetaI - cosThetaT) * Vec(n);
		return true;
	}
private:
	real nc, nt;
	Vec Fa;
	TransportMode transportMode;
};


NAMESPACE_END