#pragma once

#include "linagl.h"
#include "intersection.h"
#include "light.h"
#include "sampling.h"
#include "geometry_util.h"

NAMESPACE_BEGIN

class AreaLight : public Light {
public:
	AreaLight(const std::shared_ptr<Shape>& pShape) : shape(pShape) {}
	AreaLight(const std::shared_ptr<Shape>& pShape, const Vec3 & Lemit) : shape(pShape), Lemit(Lemit) {}
	Vec3 DirectIllumination(const Intersection &isect, const std::shared_ptr<BSDF> &bsdf,
		const Vec3 &importance, Vec3 *dir, Intersection *lightPoint, const Vec2 &u) const override {
		real pdf;
		//lightPoint->hit = shape->Sample(isect, &pdf, u);
		//lightPoint->n = lightPoint->nl = shape->GetNorm(lightPoint->hit);
		*lightPoint = shape->Sample(isect, &pdf, u);
		*dir = (lightPoint->hit - isect.hit).Norm();
		Vec3 f = bsdf->f(isect.wo, *dir);
		return importance * f * std::abs((*dir).Dot(isect.n)) * Emission() / pdf;
	}

	Vec3 Sample_Li(const Intersection &isect, Vec3 *wi, real *pdf, Intersection *lightPoint, const Vec2 &u) const override {
		//lightPoint->hit = shape->Sample(isect, pdf, u);
		//lightPoint->n = lightPoint->nl = shape->GetNorm(lightPoint->hit);
		*lightPoint = shape->Sample(isect, pdf, u);
		*wi = (lightPoint->hit - isect.hit).Norm();
		return Emission();
	}

	real Pdf_Li(const Intersection &isect, const Vec3 &wi) const override {
		return shape->Pdf(isect, wi);
	}

	Vec3 Emission() const override {
		return Lemit;
	}

	int GetId() const override {
		return shape->GetId();
	}

	Vec3 Power() const override {
		return Emission() * shape->Area() * PI;
	}

	bool IsAreaLight() const override { return true; }

	std::shared_ptr<Shape> GetShapePtr() const override { return shape; }
protected:
	//void SampleOnLight(Vec3 *pos, Vec3 *dir, Vec3 *lightNorm, real *pdfPos, real *pdfDir, const Vec2 &u, const Vec2 &v) const override{
	//	//sample a position
	//	*pos = shape->Sample(pdfPos, u);
	//	*lightNorm = shape->GetNorm(*pos);
	//	Vec3 ss, ts;
	//	CoordinateSystem(*lightNorm, &ss, &ts);
	//	Vec3 dirLocal = CosineSampleHemisphere(v);
	//	real cosTheta = dirLocal.z;
	//	*dir = (ss * dirLocal.x + ts * dirLocal.y + *lightNorm * dirLocal.z).Norm();
	//	*pdfDir = CosineHemispherePdf(cosTheta);
	//}

	void SampleOnLight(Intersection *isect, Vec3 *dir, real *pdfPos, real *pdfDir, const Vec2 &u, const Vec2 &v) const override {
		//sample a position
		*isect = shape->Sample(pdfPos, u);
		Vec3 ss, ts;
		CoordinateSystem(isect->n, &ss, &ts);
		Vec3 dirLocal = CosineSampleHemisphere(v);
		real cosTheta = dirLocal.z;
		*dir = (ss * dirLocal.x + ts * dirLocal.y + isect->n * dirLocal.z).Norm();
		*pdfDir = CosineHemispherePdf(cosTheta);
	}
private:
	std::shared_ptr<Shape> shape;
	Vec3 Lemit;
};

NAMESPACE_END