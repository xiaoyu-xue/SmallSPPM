#pragma once

#include "math/Linagl.h"
#include "common/Core.h"
#include "core/Intersection.h"
#include "math/Transform.h"
#include "core/Sampling.h"
#include "core/Microfacet.h"
#include "system/Threading.h"
#include "math/MathUtils.h"
#include "math/Frame.h"
#include "common/Core.h"

GYT_NAMESPACE_BEGIN

namespace BSDFCoordinate {
	GYT_FORCE_INLINE real CosTheta(const Vec3& w) { return w.z; }

	GYT_FORCE_INLINE real Cos2Theta(const Vec3& w) { return w.z * w.z; }

	GYT_FORCE_INLINE real AbsCosTheta(const Vec3& w) { return std::abs(w.z); }

	GYT_FORCE_INLINE real Sin2Theta(const Vec3& w) {
		return std::max((real)0, (real)1 - Cos2Theta(w));
	}

	GYT_FORCE_INLINE real SinTheta(const Vec3& w) { return std::sqrt(Sin2Theta(w)); }

	GYT_FORCE_INLINE real TanTheta(const Vec3& w) { return SinTheta(w) / CosTheta(w); }

	GYT_FORCE_INLINE real Tan2Theta(const Vec3& w) {
		return Sin2Theta(w) / Cos2Theta(w);
	}

	GYT_FORCE_INLINE real CosPhi(const Vec3& w) {
		real sinTheta = SinTheta(w);
		return (sinTheta == 0) ? 1 : Clamp(w.x / sinTheta, -1, 1);
	}

	GYT_FORCE_INLINE real SinPhi(const Vec3& w) {
		real sinTheta = SinTheta(w);
		return (sinTheta == 0) ? 0 : Clamp(w.y / sinTheta, -1, 1);
	}

	GYT_FORCE_INLINE real Cos2Phi(const Vec3& w) { return CosPhi(w) * CosPhi(w); }

	GYT_FORCE_INLINE real Sin2Phi(const Vec3& w) { return SinPhi(w) * SinPhi(w); }

	GYT_FORCE_INLINE real CosDPhi(const Vec3& wa, const Vec3& wb) {
		return Clamp(
			(wa.x * wb.x + wa.y * wb.y) / std::sqrt((wa.x * wa.x + wa.y * wa.y) *
			(wb.x * wb.x + wb.y * wb.y)),
			-1, 1);
	}
}

enum ScatterEventType
{
	BSDF_REFLECTION = 1 << 0,
	BSDF_TRANSMISSION = 1 << 1,
	BSDF_DIFFUSE = 1 << 2,
	BSDF_GLOSSY = 1 << 3,
	BSDF_SPECULAR = 1 << 4,
	BSDF_ALL_TYPES = BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_SPECULAR,
	BSDF_ALL_REFLECTION = BSDF_REFLECTION | BSDF_ALL_TYPES,
	BSDF_ALL_TRANSMISSION = BSDF_TRANSMISSION | BSDF_ALL_TYPES,
	BSDF_ALL = BSDF_ALL_REFLECTION | BSDF_ALL_TRANSMISSION
};

real FrDielectric(real cosThetaI, real etaI, real etaT);

// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
Vec3 FrConductor(real cosThetaI, const Vec3& etai, const Vec3& etat, const Vec3& k);


class Fresnel {
public:
	// Fresnel Interface
	virtual ~Fresnel(){}
	virtual Vec3 Evaluate(real cosI) const = 0;
	virtual real EvaluateScalar(real cosI) const = 0;
};


class FresnelConductor : public Fresnel {
private:
	Vec3 m_etaI, m_etaT, m_k;
public:
	FresnelConductor(const Vec3& etaI, const Vec3& etaT, const Vec3& k)
		: m_etaI(etaI), m_etaT(etaT), m_k(k) {}

	Vec3 Evaluate(real cosThetaI) const override{
		return FrConductor(std::abs(cosThetaI), m_etaI, m_etaT, m_k);
	}

	real EvaluateScalar(real cosI) const override {
		return 0;
	}
};


class FresnelDielectric : public Fresnel {
private:
	real m_etaI, m_etaT;
public:
	FresnelDielectric(real etaI, real etaT) : m_etaI(etaI), m_etaT(etaT) {}

	Vec3 Evaluate(real cosThetaI) const {
		return FrDielectric(cosThetaI, m_etaI, m_etaT);
	}

	real EvaluateScalar(real cosThetaI) const override {
		return FrDielectric(cosThetaI, m_etaI, m_etaT);
	}
};


class BxDF;
class BSDF {
protected:
	Vec3 mGeometryNormal, mNormal;
	Frame mFrame;
	static constexpr int MaxBSDFs = 8;
	BxDF* mBxdfs[MaxBSDFs];
	int mBsdfCnt = 0;
	bool mIsDelta = false;
public:
	BSDF(const Intersection &isect)
		: mNormal(isect.mNormal), mGeometryNormal(isect.mGeometryNormal), mFrame(isect.mShadingDpDu, isect.mShadingDpDv, isect.mNormal)
	{
	}
	GYT_FORCE_INLINE Vec3 WorldToLocal(const Vec3& v) const 
	{
		return mFrame.WorldToLocal(v);
	}
	GYT_FORCE_INLINE Vec3 LocalToWorld(const Vec3& v) const
	{
		return mFrame.LocalToWorld(v);
	}
	~BSDF(){}
	void Add(BxDF* bsdf);
	real Pdf(const Vec3 &wo, const Vec3 &wi, ScatterEventType flags = BSDF_ALL) const;
	Vec3 Sample(const Vec3 &wo, Vec3 *wi, real *pdf, const Vec3 &rand = Vec3(0, 0, 0), ScatterEventType flags = BSDF_ALL) const;
	Vec3 Evaluate(const Vec3 &wo, const Vec3 &wi, ScatterEventType flags = BSDF_ALL) const;
	bool IsDelta() const { 
		return mIsDelta;
	}
	bool MatchScatterType(ScatterEventType flag) const;

};


class BxDF {
public:
	BxDF(ScatterEventType scatterType) : scatterEventType(scatterType){}
	virtual ~BxDF(){}
	virtual real Pdf(const Vec3& wo, const Vec3& wi) const = 0;
	virtual Vec3 Sample(const Vec3& wo, Vec3* wi, real* pdf, const Vec3& rand = Vec3(0, 0, 0)) const = 0;
	virtual Vec3 Evaluate(const Vec3& wo, const Vec3& wi) const = 0;
	virtual bool IsDelta() const { return scatterEventType & BSDF_SPECULAR; }
	virtual bool MatchesTypes(ScatterEventType flags) const { return (scatterEventType & flags) == scatterEventType; }
	const ScatterEventType scatterEventType;
};

class DiffuseBSDF : public BxDF {
private:
	Vec3 mR;
public:
	DiffuseBSDF(Vec3 r) 
		: BxDF(ScatterEventType(BSDF_REFLECTION | BSDF_DIFFUSE)), mR(r) 
	{ 
	}

	real Pdf(const Vec3& wo, const Vec3& wi) const override 
	{
		return SameHemisphere(wo, wi) ? BSDFCoordinate::AbsCosTheta(wi) * INV_PI : 0.f;
	}

	Vec3 Sample(const Vec3& wo, Vec3* wi, real* pdf, const Vec3& rand) const override 
	{
		*wi = CosineSampleHemisphere(Vec2(rand[0], rand[1]));
		if (wo.z < 0) wi->z *= -1;
		*pdf = Pdf(wo, *wi);
		return Evaluate(wo, *wi);
	}

	Vec3 Evaluate(const Vec3 &wo, const Vec3 &wi) const override {
		return mR * INV_PI;
	}
};


class SpecularBSDF : public BxDF {
public:
	SpecularBSDF(Vec3 r = Vec3(1.0, 1.0, 1.0)) 
		: BxDF(ScatterEventType(BSDF_ALL_REFLECTION | BSDF_SPECULAR)), R(r) 
	{
	}

	real Pdf(const Vec3 &wo, const Vec3 &wi) const override 
	{
		return 0.0;
	}

	Vec3 Sample(const Vec3& wo, Vec3* wi, real* pdf, const Vec3& rand) const override 
	{
		*wi = Vec3(-wo.x, -wo.y, wo.z);
		*pdf = 1.0;
		real cosTheta = BSDFCoordinate::AbsCosTheta(*wi);
		return R / cosTheta;
	}

	Vec3 Evaluate(const Vec3 &wo, const Vec3 &wi) const override 
	{
		return Vec3();
	}

	inline bool IsDelta() const override { return true; }
private:
	Vec3 R;
};

//TODO

class TransmissionBSDF : public BxDF 
{
private:
	real m_nc, m_nt;
	Vec3 mR, mT;
	TransportMode mTransportMode;
public:
	TransmissionBSDF(const Vec3 kt = Vec3(1.0, 1.0, 1.0), const Vec3 kr = Vec3(1.0, 1.0, 1.0), 
		TransportMode mode = TransportMode::Radiance, real eta1 = 1.0, real eta2 = 1.5) 
		:BxDF(ScatterEventType(BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_SPECULAR)),
		mT(kt), mR(kr), m_nc(eta1), m_nt(eta2), mTransportMode(mode) 
	{
	}

	real Pdf(const Vec3 &wo, const Vec3 &wi) const override 
	{
		return 0.0;
	}

	Vec3 Sample(const Vec3& wo, Vec3* wi, real* pdf, const Vec3& rand) const override 
	{

		real Re = FrDielectric(BSDFCoordinate::CosTheta(wo), m_nc, m_nt);

		bool entering = BSDFCoordinate::CosTheta(wo) > 0;
		real etaI = entering ? m_nc : m_nt;
		real etaT = entering ? m_nt : m_nc;
		Vec3 nForward = BSDFCoordinate::CosTheta(wo) > 0 ? Vec3(0, 0, 1) : Vec3(0, 0, -1);
		real P = Re * 0.5f + 0.25f;
		P = Re;
		//P = 1.0;
		if (rand[1] < P) {
			//*wi = (2 * wo.Dot(nForward) * nForward - wo).Norm();
			*wi = Vec3(-wo.x, -wo.y, wo.z);
			*pdf = P;
			real cosTheta = BSDFCoordinate::AbsCosTheta(*wi);
			return mR * Re / cosTheta;
		}
		else {
			if (!Refract(wo, nForward, etaI / etaT, wi)) {
				return Vec3();
			}
			*pdf = 1 - P;
			real cosTheta = BSDFCoordinate::AbsCosTheta(*wi);
			real eta = etaI / etaT;
			if (mTransportMode == TransportMode::Radiance) {
				return eta * eta * mT * (1.f - Re) / cosTheta;
			}
			else {
				return mT * (1.f - Re) / cosTheta;
			}
		}
	}

	Vec3 Evaluate(const Vec3 &wo, const Vec3 &wi) const override 
	{
		return Vec3();
	}

	real Fresnell(const Vec3 &wo, const Vec3 &td, const Vec3 &n, const Vec3 &nl) const 
	{
		bool into = (n.Dot(nl) > 0.0);
		real nnt = into ? m_nc / m_nt : m_nt / m_nc, ddn = (-1 * wo).Dot(nl);//, cos2t;
		//cos2t = 1 - nnt * nnt*(1 - ddn * ddn);
		real a = m_nt - m_nc, b = m_nt + m_nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : td.Dot(n));
		real Re = R0 + (1 - R0) * c * c* c * c * c;
		return Re;
	}

	inline bool IsDelta() const override { return true; }

	static bool Refract(const Vec3 &wi, const Vec3 &n, real eta, Vec3 *wt) 
	{
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
};


class MicrofacetReflectionBSDF : public BxDF {
private:
	std::unique_ptr<MicrofacetDistribution> mpDistribution;
	std::unique_ptr<Fresnel> mpFresnel;
	Vec3 mR;

public:

	MicrofacetReflectionBSDF(MicrofacetDistribution* distribution, Fresnel *fresnel, const Vec3& R) 
		:BxDF(ScatterEventType(BSDF_REFLECTION | BSDF_GLOSSY)), mpDistribution(distribution), mpFresnel(fresnel), mR(R) 
	{
	}

	real Pdf(const Vec3& wo, const Vec3& wi) const override 
	{
		real pdf = 0;
		Vec3 wh;
		real dwh_dwi;
		wh = (wo + wi).Norm();
		dwh_dwi = 1.f / (4 * Dot(wi, wh));
		wh = wh * Sgn(BSDFCoordinate::CosTheta(wh));
		pdf = std::abs(mpDistribution->Pdf(wo, wh) * dwh_dwi);
		return pdf;
	}

	Vec3 Sample(const Vec3& wo, Vec3* wi, real* pdf, const Vec3& rand) const override 
	{
		if (wo.z == 0) return 0.;
		//Vec3 wh = distribution->Sample_wh(wo, Vec2(rand[0], rand[1]));
		Vec3 wh = mpDistribution->Sample_wh(Sgn(BSDFCoordinate::CosTheta(wo)) * wo, Vec2(rand[0], rand[1]));
		Vec3 microfacetPdf = mpDistribution->Pdf(wo, wh);;
		if (microfacetPdf == 0) return Vec3();

		//if (Dot(wo, wh) < 0) return 0.;   // Should be rare
		*wi = Reflect(wo, wh);
		if (!SameHemisphere(wo, *wi)) return Vec3(0, 0, 0);

		wh = Sgn(BSDFCoordinate::CosTheta(wh)) * wh;

		// Compute PDF of _wi_ for microfacet reflection
		*pdf = std::abs(mpDistribution->Pdf(wo, wh) / (4 * Dot(*wi, wh)));
		return Evaluate(wo, *wi);
	}

	Vec3 Evaluate(const Vec3& wo, const Vec3& wi) const override {
		if (!SameHemisphere(wo, wi)) return Vec3();

		real cosThetaO = BSDFCoordinate::AbsCosTheta(wo), cosThetaI = BSDFCoordinate::AbsCosTheta(wi);
		Vec3 wh = wi + wo;
		// Handle degenerate cases for microfacet reflection
		if (cosThetaI == 0 || cosThetaO == 0) return Vec3(0, 0, 0);
		if (wh.x == 0 && wh.y == 0 && wh.z == 0) return Vec3(0, 0, 0);
		wh = wh.Norm();
		// For the Fresnel call, make sure that wh is in the same hemisphere
		// as the surface normal, so that TIR is handled correctly.
		wh = Sgn(BSDFCoordinate::CosTheta(wh)) * wh;
		Vec3 F = mpFresnel->Evaluate(Dot(wo, wh));
		Vec3 D = mpDistribution->D(wh);
		Vec3 G = mpDistribution->G(wo, wi, wh);
		return mR * F * D * G / (4 * cosThetaI * cosThetaO);

	}
};


class MicrofacetTransmissionBSDF : public BxDF {
private:
	const Vec3 mT;
	const real m_etaA, m_etaB;
	const FresnelDielectric mFresnel;
	const TransportMode mTransportMode;
	std::unique_ptr<MicrofacetDistribution> mpDistribution;
public:
	// MicrofacetTransmission Public Methods
	MicrofacetTransmissionBSDF(MicrofacetDistribution* distribution, const Vec3& T, real etaA, real etaB, TransportMode mode)
		: BxDF(ScatterEventType(BSDF_TRANSMISSION | BSDF_GLOSSY)),
		mT(T),
		mpDistribution(distribution),
		m_etaA(etaA),
		m_etaB(etaB),
		mFresnel(etaA, etaB),
		mTransportMode(mode)
	{
	}

	Vec3 Evaluate(const Vec3& wo, const Vec3& wi) const override 
	{
		if (SameHemisphere(wo, wi)) return Vec3();  // transmission only

		real cosThetaO = BSDFCoordinate::CosTheta(wo);
		real cosThetaI = BSDFCoordinate::CosTheta(wi);
		if (cosThetaI == 0 || cosThetaO == 0) return Vec3();

		real eta = BSDFCoordinate::CosTheta(wo) > 0 ? (m_etaB / m_etaA) : (m_etaA / m_etaB);
		Vec3 wh = (wo + wi * eta).Norm();
		if (wh.z < 0) wh = -wh;

		Vec3 F = mFresnel.Evaluate(Dot(wo, wh));

		real sqrtDenom = Dot(wo, wh) + eta * Dot(wi, wh);
		real factor = (mTransportMode == TransportMode::Radiance) ? (1 / eta) : 1;

		return (Vec3(1.f) - F) * mT *
			std::abs(mpDistribution->D(wh) * mpDistribution->G(wo, wi, wh) * eta * eta *
				std::abs(Dot(wi, wh)) * std::abs(Dot(wo, wh)) * factor * factor /
				(cosThetaI * cosThetaO * sqrtDenom * sqrtDenom));
	}

	Vec3 Sample(const Vec3& wo, Vec3* wi, real* pdf, const Vec3& rand) const override 
	{
		if (wo.z == 0) return 0.;
		//Vec3 wh = distribution->Sample_wh(wo, Vec2(rand[0], rand[1]));
		Vec3 wh = mpDistribution->Sample_wh(Sgn(BSDFCoordinate::CosTheta(wo)) * wo, Vec2(rand[0], rand[1]));
		real eta = BSDFCoordinate::CosTheta(wo) > 0 ? (m_etaA / m_etaB) : (m_etaB / m_etaA);
		if (!TransmissionBSDF::Refract(wo, Faceforward(wh, wo), eta, wi)) {
			*pdf = 0;
			return Vec3();
		}
		*pdf = Pdf(wo, *wi);

		return Evaluate(wo, *wi);
	}

	real Pdf(const Vec3& wo, const Vec3& wi) const override
	{
		if (SameHemisphere(wo, wi)) return 0;

		real eta = BSDFCoordinate::CosTheta(wo) > 0 ? (m_etaB / m_etaA) : (m_etaA / m_etaB);
		Vec3 wh = -((wo + wi * eta).Norm());

		wh = Sgn(BSDFCoordinate::CosTheta(wh)) * wh;
		// Compute change of variables _dwh\_dwi_ for microfacet transmission
		real sqrtDenom = Dot(wo, wh) + eta * Dot(wi, wh);
		real dwh_dwi = std::abs((eta * eta * Dot(wi, wh)) / (sqrtDenom * sqrtDenom));

		return mpDistribution->Pdf(wo, wh) * dwh_dwi;
	}
};


class RoughDielectricBSDF : public BxDF{
private:
	const Vec3 mR, mT;
	const real m_etaA, m_etaB;
	const FresnelDielectric mFresnel;
	const TransportMode mTransportMode;
	std::unique_ptr<MicrofacetDistribution> distribution;

public:
	RoughDielectricBSDF(MicrofacetDistribution* distribution, const Vec3 &R, const Vec3& T, real etaA, real etaB, TransportMode mode)
		: BxDF(ScatterEventType(BSDF_REFLECTION | BSDF_TRANSMISSION | BSDF_GLOSSY)),
		mR(R), mT(T),
		distribution(distribution),
		m_etaA(etaA),
		m_etaB(etaB),
		mFresnel(etaA, etaB),
		mTransportMode(mode) 
	{
	}

	Vec3 Evaluate(const Vec3& wo, const Vec3& wi) const override 
	{
		bool reflect = BSDFCoordinate::CosTheta(wo) * BSDFCoordinate::CosTheta(wi) > 0;
		Vec3 wh;
		if (reflect) 
		{
			real cosThetaO = BSDFCoordinate::AbsCosTheta(wo), cosThetaI = BSDFCoordinate::AbsCosTheta(wi);
			wh = wi + wo;
			// Handle degenerate cases for microfacet reflection
			if (cosThetaI == 0 || cosThetaO == 0) return Vec3(0, 0, 0);
			if (wh.x == 0 && wh.y == 0 && wh.z == 0) return Vec3(0, 0, 0);
			wh = wh.Norm();
			// For the Fresnel call, make sure that wh is in the same hemisphere
			// as the surface normal, so that TIR is handled correctly.
			wh = Sgn(BSDFCoordinate::CosTheta(wh)) * wh;
			real F = mFresnel.EvaluateScalar(Dot(wo, wh));
			Vec3 D = distribution->D(wh);
			Vec3 G = distribution->G(wo, wi, wh);
			return mR * F * D * G / (4 * cosThetaI * cosThetaO);
		}
		else 
		{

			real cosThetaO = BSDFCoordinate::CosTheta(wo);
			real cosThetaI = BSDFCoordinate::CosTheta(wi);
			if (cosThetaI == 0 || cosThetaO == 0) return Vec3();

			real eta = BSDFCoordinate::CosTheta(wo) > 0 ? (m_etaB / m_etaA) : (m_etaA / m_etaB);
			Vec3 wh = (-(wo + wi * eta).Norm());
			wh = Sgn(BSDFCoordinate::CosTheta(wh)) * wh;

			real F = mFresnel.EvaluateScalar(Dot(wo, wh));
			real D = distribution->D(wh);
			real G = distribution->G(wo, wi, wh);


			real sqrtDenom = Dot(wo, wh) + eta * Dot(wi, wh);
			real factor = (mTransportMode == TransportMode::Radiance) ? (1 / eta) : 1;
			
			real value = (1 - F) * D * G * eta * eta * std::abs(Dot(wi, wh)) * std::abs(Dot(wo, wh)) / std::abs(cosThetaI * cosThetaO * sqrtDenom * sqrtDenom);

			return mT * std::abs(value * factor * factor);

		}
	}

	real Pdf(const Vec3& wo, const Vec3& wi) const override 
	{

		bool reflect = BSDFCoordinate::CosTheta(wo) * BSDFCoordinate::CosTheta(wi) > 0;
		real pdf = 0;
		Vec3 wh;
		real dwh_dwi ;
		if (reflect) {
			wh = (wo + wi).Norm();
			dwh_dwi = 1.f / (4 * Dot(wi, wh));
		}
		else {
			real eta = BSDFCoordinate::CosTheta(wo) > 0 ? (m_etaB / m_etaA) : (m_etaA / m_etaB);
			wh = -((wo + wi * eta).Norm());

			// Compute change of variables _dwh\_dwi_ for microfacet transmission
			real sqrtDenom = Dot(wo, wh) + eta * Dot(wi, wh);
			dwh_dwi = std::abs((eta * eta * Dot(wi, wh)) / (sqrtDenom * sqrtDenom));
		}

		wh = wh * Sgn(BSDFCoordinate::CosTheta(wh));

		pdf = std::abs(distribution->Pdf(wo, wh) * dwh_dwi);

		real F = mFresnel.EvaluateScalar(Dot(wo, wh));

		pdf *= reflect ? F : (1 - F);

		return pdf;
	}

	Vec3 Sample(const Vec3& wo, Vec3* wi, real* pdf, const Vec3& rand) const override 
	{
		Vec3 wh = distribution->Sample_wh(Sgn(BSDFCoordinate::CosTheta(wo)) * wo, Vec2(rand[0], rand[1]));
		real microfacetPdf = distribution->Pdf(wo, wh);
		if (microfacetPdf == 0) return Vec3();

		real F = mFresnel.EvaluateScalar(Dot(wo, wh));
		real P = F;
		*pdf = 0;
		DEBUG_PIXEL_IF(ThreadIndex()) {
			std::cout << "Sample_f fresnel : " << F << std::endl;
		}
		if (rand[2] < P) {
			if (wo.z == 0) return Vec3();
			//if (Dot(wo, wh) < 0) return 0.;   // Should be rare
			*wi = Reflect(wo, wh);

			if (!SameHemisphere(wo, *wi)) return Vec3(0, 0, 0);
			
			real dwh_dwi = 1.f / (4 * Dot(*wi, wh));

			*pdf = P * std::abs(distribution->Pdf(wo, wh) * dwh_dwi);

			real cosThetaO = BSDFCoordinate::AbsCosTheta(wo), cosThetaI = BSDFCoordinate::AbsCosTheta(*wi);

			Vec3 D = distribution->D(wh);

			Vec3 G = distribution->G(wo, *wi, wh);

			return mR * F * D * G / (4 * cosThetaI * cosThetaO);
		}
		else 
		{

			DEBUG_PIXEL_IF(ThreadIndex()) {
				std::cout << "Sample_f wo : " << wo << " Dot(wo, wh): " << Dot(wo, wh) << std::endl;
			}

			if (wo.z == 0) return Vec3();
			//if (Dot(wo, wh) < 0) return 0.;  // Should be rare

			real etaRefract = BSDFCoordinate::CosTheta(wo) > 0 ? (m_etaA / m_etaB) : (m_etaB / m_etaA);
			if (!TransmissionBSDF::Refract(wo, Faceforward(wh, wo), etaRefract, wi)) {
				*pdf = 0;
				return Vec3();
			}

			if (BSDFCoordinate::CosTheta(wo) * BSDFCoordinate::CosTheta(*wi) >= 0) return Vec3();

			real eta = BSDFCoordinate::CosTheta(wo) > 0 ? (m_etaB / m_etaA) : (m_etaA / m_etaB);
			real sqrtDenom = Dot(wo, wh) + eta * Dot(*wi, wh);
			real dwh_dwi = std::abs((eta * eta * Dot(*wi, wh)) / (sqrtDenom * sqrtDenom));
			*pdf = std::abs((1 - P) * distribution->Pdf(wo, wh) * dwh_dwi);


			real cosThetaO = BSDFCoordinate::CosTheta(wo);
			real cosThetaI = BSDFCoordinate::CosTheta(*wi);
			if (cosThetaI == 0 || cosThetaO == 0) return Vec3();

			real D = distribution->D(wh);

			real G = distribution->G(wo, *wi, wh);

			real factor = (mTransportMode == TransportMode::Radiance) ? (1 / eta) : 1;

			real value = (1 - F) * D * G * eta * eta * std::abs(Dot(*wi, wh)) * std::abs(Dot(wo, wh)) / std::abs(cosThetaI * cosThetaO * sqrtDenom * sqrtDenom);

			return mT * std::abs(value * factor * factor);
		}

	}
};

GYT_NAMESPACE_END