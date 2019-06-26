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
	Vec DirectIllumination(const Intersection &isect, const std::shared_ptr<BSDF> &bsdf,
		const Vec &importance, Vec *dir, Intersection *lightPoint, Vec u) const {
		real pdf;
		lightPoint->hit = shape->Sample(isect, &pdf, u);
		lightPoint->n = lightPoint->nl = shape->GetNorm(lightPoint->hit);
		*dir = (lightPoint->hit - isect.hit).norm();
		Vec f = bsdf->f(isect.wo, *dir);
		return importance * f * std::abs((*dir).dot(isect.n)) * Emission() / pdf;
	}

	Vec Sample_Li(const Intersection &isect, Vec *wi, real *pdf, Intersection *lightPoint, Vec u) const {
		lightPoint->hit = shape->Sample(isect, pdf, u);
		lightPoint->n = lightPoint->nl = shape->GetNorm(lightPoint->hit);
		*wi = (lightPoint->hit - isect.hit).norm();
		return Emission();
	}

	real Pdf_Li(const Intersection &isect, const Vec &wi) const {
		return shape->Pdf(isect, wi);
	}

	Vec Emission() const {
		return shape->GetEmission();
	}

	int GetId() const {
		return shape->GetId();
	}

	Vec Power() const {
		return Emission() * shape->Area() * PI;
	}

	bool IsAreaLight() const { return true; }

	std::shared_ptr<Shape> GetShapePtr() const { return shape; }
protected:
	void SampleOnLight(Vec *pos, Vec *dir, Vec *lightNorm, real *pdfPos, real *pdfDir, Vec u, Vec v) const {
		//sample a position
		*pos = shape->Sample(pdfPos, u);
		*lightNorm = shape->GetNorm(*pos);
		Vec ss, ts;
		CoordinateSystem(*lightNorm, &ss, &ts);
		Vec dirLocal = CosineSampleHemisphere(v);
		real cosTheta = dirLocal.z;
		*dir = (ss * dirLocal.x + ts * dirLocal.y + *lightNorm * dirLocal.z).norm();
		*pdfDir = CosineHemispherePdf(cosTheta);
	}
private:
	std::shared_ptr<Shape> shape;
};

NAMESPACE_END