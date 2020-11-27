#pragma once
#include "Light.h"

NAMESPACE_BEGIN

class PointLight : public Light {
public:
	PointLight(const Vec3& p, const Vec3& intensity) : mLightPosition(p), mI(intensity) {}

	Vec3 Emission() const override { return mI; }

	bool IsDeltaLight() const override { return true; }

	Vec3 Sample_Li(const Intersection& isect, Vec3* wi, real* pdf, Intersection* lightPoint, const Vec2& u) const override;

	real Pdf_Li(const Intersection& isect, const Vec3& wi) const override 
	{
		return 0;
	}

	FORCE_INLINE Vec3 Power() const override 
	{
		return 4 * PI * mI;
	}

	Vec3 SampleLight(Intersection* isect, Vec3* dir, real* pdfPos, real* pdfDir, const Vec2& u, const Vec2& v) const override;

	Vec3 SampleOnePoint(Intersection* isect, real* pdf, const Vec2& u) const override 
	{
		isect->hit = mLightPosition;
		*pdf = 1.f;
		return mI;
	}
private:
	Vec3 mLightPosition;
	Vec3 mI;
};

NAMESPACE_END