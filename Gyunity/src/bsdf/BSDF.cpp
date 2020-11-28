#include "BSDF.h"

GY_NAMESPACE_BEGIN

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
	bxdfs[nBSDFs++] = bsdf;
	isDelta |= (bsdf->scatterEventType & BSDF_SPECULAR);
}

real BSDF::Pdf(const Vec3& wo, const Vec3& wi, ScatterEventType flag) const {
	if (nBSDFs == 0) return 0.f;
	Vec3 woLocal = WorldToLocal(wo);
	Vec3 wiLocal = WorldToLocal(wi);
	real pdf = 0.f;
	int matchingComps = 0;
	for (int i = 0; i < nBSDFs; ++i) {
		if (bxdfs[i]->MatchesTypes(flag)) {
			++matchingComps;
			pdf += bxdfs[i]->Pdf(woLocal, wiLocal);
		}
	}
	//DEBUG_PIXEL_IF(ThreadIndex()) {
	//	std::cout << "Pdf: " << pdf << std::endl;
	//}
	return matchingComps > 0 ? pdf / matchingComps : 0.f;
}

Vec3 BSDF::Sample_f(const Vec3& wo, Vec3* wi, real* pdf, const Vec3& rand, ScatterEventType flags) const {
	int sampleIndex = rand[0] * nBSDFs;
	Vec3 woLocal = WorldToLocal(wo);
	Vec3 wiLocal;
	BxDF* bxdf = bxdfs[sampleIndex];

	*pdf = 0;

	//DEBUG_PIXEL_IF(ThreadIndex()) {
	//	std::cout << "woLocal: " << woLocal << std::endl;
	//}

	Vec3 f = bxdf->Sample_f(woLocal, &wiLocal, pdf, Vec3(rand[0], rand[1], rand[2]));
	*wi = LocalToWorld(wiLocal);

	//DEBUG_PIXEL_IF(ThreadIndex()) {
	//	std::cout << "BSDF(0): " << f << " Pdf: " << *pdf << std::endl;
	//}

	for (int i = 0; i < nBSDFs; ++i) {
		if (i != sampleIndex) {
			*pdf += bxdfs[i]->Pdf(woLocal, wiLocal);
		}
	}
	*pdf /= nBSDFs;

	if (!bxdf->IsDelta()) {
		bool reflect = Dot(*wi, ng) * Dot(wo, ng) > 0;
		//f = Vec3(0, 0, 0);
		//DEBUG_PIXEL_IF(ThreadIndex()) {
		//	std::cout << "BSDF(1): " << f << " Pdf: " << *pdf << std::endl;
		//}
		for (int i = 0; i < nBSDFs; ++i) {
			if((reflect && (bxdfs[i]->scatterEventType & BSDF_REFLECTION)) ||
				(!reflect && (bxdfs[i]->scatterEventType & BSDF_TRANSMISSION)))
				if(i != sampleIndex) f += bxdfs[i]->f(woLocal, wiLocal);
				//f += bxdfs[i]->f(woLocal, wiLocal);
		}
	}
	//DEBUG_PIXEL_IF(ThreadIndex()) {
	//	std::cout << "BSDF(2): " << f << " Pdf: " << *pdf << std::endl;
	//	if (std::isnan((*wi).x)) {
	//		std::cout << "wi: " << *wi << " wo: " << wo << std::endl;
	//	}
	//}
	return f;
}

Vec3 BSDF::f(const Vec3& wo, const Vec3& wi, ScatterEventType flags) const{
	Vec3 f(0, 0, 0);
	Vec3 woLocal = WorldToLocal(wo);
	Vec3 wiLocal = WorldToLocal(wi);
	bool reflect = Dot(wi, ng) * Dot(wo, ng) > 0;
	//DEBUG_PIXEL_IF(ThreadIndex()) {
	//	std::cout << "Reflect?: " << reflect << " " << Dot(wi, ng)  << " " << Dot(wo, ng) << std::endl;
	//}
	for (int i = 0; i < nBSDFs; ++i) {
		if ((reflect && (bxdfs[i]->scatterEventType & BSDF_REFLECTION)) ||
			(!reflect && (bxdfs[i]->scatterEventType & BSDF_TRANSMISSION)))
			f += bxdfs[i]->f(woLocal, wiLocal);
	}
	//DEBUG_PIXEL_IF(ThreadIndex()) {
	//	std::cout << "BSDFs: " << f << std::endl;
	//}
	return f;
}


bool BSDF::MatchScatterType(ScatterEventType flag) const {
	int ret = 0;
	for (int i = 0; i < nBSDFs; ++i) {
		ret |= (bxdfs[i]->scatterEventType & flag);
	}
	return ret > 0;
}
GY_NAMESPACE_END