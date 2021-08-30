#pragma once

#include "math/Linagl.h"
#include "visual/Intersection.h"
#include "light/Light.h"
#include "visual/Sampling.h"
#include "math/GeometryUtils.h"
#include "shape/Shape.h"

GYT_NAMESPACE_BEGIN

class AreaLight : public Light 
{
private:
	std::shared_ptr<Shape> mpShape;
	Vec3 mLemit;
public:
	AreaLight(const std::shared_ptr<Shape>& pShape) : mpShape(pShape) {}
	AreaLight(const std::shared_ptr<Shape>& pShape, const Vec3 & Lemit) : mpShape(pShape), mLemit(Lemit) {}

	Vec3 Sample_Li(const Intersection &isect, Vec3 *wi, real *pdf, Intersection *lightPoint, const Vec2 &u) const override 
	{
		*lightPoint = mpShape->Sample(isect, pdf, u);
		*wi = (lightPoint->mPos - isect.mPos).Norm();
		return Emission(*lightPoint, -*wi);
	}

	real Pdf_Li(const Intersection &isect, const Vec3 &wi) const override 
	{
		return mpShape->Pdf(isect, wi);
	}

	Vec3 Emission() const override 
	{
		return mLemit;
	}

	Vec3 Emission(const Intersection &isect, const Vec3 &w) const override 
	{
		if (Dot(isect.mNormal, w) > 0) return mLemit;
		else return Vec3(0, 0, 0);

	}

	Vec3 Power() const override 
	{
		return Emission() * mpShape->Area() * PI;
	}

	bool IsAreaLight() const override { return true; }

	std::shared_ptr<Shape> GetShapePtr() const override { return mpShape; }

	Vec3 SampleOnePoint(Intersection* isect, real* pdf, const Vec2& u) const override 
	{
		*isect = mpShape->Sample(pdf, u);
		return mLemit;
	}

	void SampleOnLight(Intersection* isect, Vec3* dir, real* pdfPos, real* pdfDir, const Vec2& u, const Vec2& v) const override
	{
		//sample a position
		*isect = mpShape->Sample(pdfPos, u);
		Vec3 ss, ts;
		CoordinateSystem(isect->mNormal, &ss, &ts);
		Vec3 dirLocal = CosineSampleHemisphere(v);
		real cosTheta = dirLocal.z;
		*dir = (ss * dirLocal.x + ts * dirLocal.y + isect->mNormal * dirLocal.z).Norm();
		*pdfDir = CosineHemispherePdf(cosTheta);
	}

protected:



};

GYT_NAMESPACE_END