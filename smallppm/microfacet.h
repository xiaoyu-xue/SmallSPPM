#pragma once

#include "def.h"
#include "linagl.h"

NAMESPACE_BEGIN

class MicrofacetDistribution {
public:
	virtual ~MicrofacetDistribution() { }
	virtual real D(const Vec3 &wh) const = 0;
	virtual real Lambda(const Vec3& w) const = 0;
	virtual real G1(const Vec3& w, const Vec3 &wh) const {
		return 1 / (1 + Lambda(w));
	}
	virtual real G(const Vec3& wo, const Vec3& wi, const Vec3 &wh) const {
		return 1 / (1 + Lambda(wo) + Lambda(wi));
	}
	virtual Vec3 Sample_wh(const Vec3& wo, const Vec2& u) const = 0;
	virtual real Pdf(const Vec3& wo, const Vec3& wh) const;

protected:
	MicrofacetDistribution(bool sampleVisibleArea) : sampleVisibleArea(sampleVisibleArea) {}

	const bool sampleVisibleArea;

};

//class BeckmannDistribution : public MicrofacetDistribution {
//public:
//	BeckmannDistribution(real alphax, real alphay, bool samplevis = true)
//		: MicrofacetDistribution(samplevis), alphax(alphax), alphay(alphay) {}
//	real D(const Vec3 &wh) const;
//	Vec3 Sample_wh(const Vec3 &wo, const Vec2 &u) const;
//private:
//	real Lambda(const Vec3& w) const;
//
//	real alphax, alphay;
//};

class GGXDistribution : public MicrofacetDistribution {
public:
	GGXDistribution(real alpha, bool samplevis = false) : MicrofacetDistribution(samplevis), alpha(alpha) {

	}
	real D(const Vec3& wh) const;
	real G1(const Vec3& v, const Vec3 &wh) const;
	real G(const Vec3& wo, const Vec3& wi, const Vec3 &wh) const;
	Vec3 Sample_wh(const Vec3& wo, const Vec2& u) const;
	real Lambda(const Vec3& w) const {
		return 0;
	}
private:
	real alpha;
};

NAMESPACE_END