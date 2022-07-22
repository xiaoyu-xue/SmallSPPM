#include "BSDF.h"

GYT_NAMESPACE_BEGIN

real FrDielectric(real cosThetaI, real etaI, real etaT) {
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

// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
Vec3 FrConductor(real cosThetaI, const Vec3& etai, const Vec3& etat, const Vec3& k) {
	cosThetaI = Clamp(cosThetaI, -1, 1);
	Vec3 eta = etat / etai;
	Vec3 etak = k / etai;

	real cosThetaI2 = cosThetaI * cosThetaI;
	real sinThetaI2 = 1. - cosThetaI2;
	Vec3 eta2 = eta * eta;
	Vec3 etak2 = etak * etak;

	Vec3 t0 = eta2 - etak2 - sinThetaI2;
	Vec3 a2plusb2 = Sqrt(t0 * t0 + 4 * eta2 * etak2);
	Vec3 t1 = a2plusb2 + cosThetaI2;
	Vec3 a = Sqrt(0.5f * (a2plusb2 + t0));
	Vec3 t2 = (real)2 * cosThetaI * a;
	Vec3 Rs = (t1 - t2) / (t1 + t2);

	Vec3 t3 = cosThetaI2 * a2plusb2 + sinThetaI2 * sinThetaI2;
	Vec3 t4 = t2 * sinThetaI2;
	Vec3 Rp = Rs * (t3 - t4) / (t3 + t4);

	return 0.5 * (Rp + Rs);
}

void BSDF::Add(BxDF* bsdf) {
	mBxdfs[mBsdfCnt++] = bsdf;
	mIsDelta |= (bsdf->scatterEventType & BSDF_SPECULAR);
}

real BSDF::Pdf(const Vec3& wo, const Vec3& wi, ScatterEventType flag) const {
	if (mBsdfCnt == 0) return 0.f;
	Vec3 woLocal = WorldToLocal(wo);
	Vec3 wiLocal = WorldToLocal(wi);
	real pdf = 0.f;
	int matchingComps = 0;
	for (int i = 0; i < mBsdfCnt; ++i) {
		if (mBxdfs[i]->MatchesTypes(flag)) {
			++matchingComps;
			pdf += mBxdfs[i]->Pdf(woLocal, wiLocal);
		}
	}

	return matchingComps > 0 ? pdf / matchingComps : 0.f;
}

Vec3 BSDF::Sample(const Vec3& wo, Vec3* wi, real* pdf, const Vec3& rand, ScatterEventType flags) const {
	int sampleIndex = rand[0] * mBsdfCnt;
	Vec3 woLocal = WorldToLocal(wo);
	Vec3 wiLocal;
	BxDF* bxdf = mBxdfs[sampleIndex];

	*pdf = 0;

	Vec3 f = bxdf->Sample(woLocal, &wiLocal, pdf, Vec3(rand.x, rand.y, rand.z));
	*wi = LocalToWorld(wiLocal);

	for (int i = 0; i < mBsdfCnt; ++i) {
		if (i != sampleIndex) {
			*pdf += mBxdfs[i]->Pdf(woLocal, wiLocal);
		}
	}
	*pdf /= mBsdfCnt;

	if (!bxdf->IsDelta()) {
		bool reflect = Dot(*wi, mGeometryNormal) * Dot(wo, mGeometryNormal) > 0;

		for (int i = 0; i < mBsdfCnt; ++i) {
			if((reflect && (mBxdfs[i]->scatterEventType & BSDF_REFLECTION)) ||
				(!reflect && (mBxdfs[i]->scatterEventType & BSDF_TRANSMISSION)))
				if(i != sampleIndex) f += mBxdfs[i]->Evaluate(woLocal, wiLocal);
				//f += bxdfs[i]->f(woLocal, wiLocal);
		}
	}

	return f;
}

Vec3 BSDF::Evaluate(const Vec3& wo, const Vec3& wi, ScatterEventType flags) const{
	Vec3 f(0, 0, 0);
	Vec3 woLocal = WorldToLocal(wo);
	Vec3 wiLocal = WorldToLocal(wi);
	bool reflect = Dot(wi, mGeometryNormal) * Dot(wo, mGeometryNormal) > 0;

	for (int i = 0; i < mBsdfCnt; ++i) {
		if ((reflect && (mBxdfs[i]->scatterEventType & BSDF_REFLECTION)) ||
			(!reflect && (mBxdfs[i]->scatterEventType & BSDF_TRANSMISSION)))
			f += mBxdfs[i]->Evaluate(woLocal, wiLocal);
	}

	return f;
}


bool BSDF::MatchScatterType(ScatterEventType flag) const {
	int ret = 0;
	for (int i = 0; i < mBsdfCnt; ++i) {
		ret |= (mBxdfs[i]->scatterEventType & flag);
	}
	return ret > 0;
}
GYT_NAMESPACE_END