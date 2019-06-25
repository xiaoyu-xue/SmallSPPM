#pragma once

#include "linagl.h"
#include "material.h"
#include "ray.h"
#include "intersection.h"
#include "bsdf.h"

class Shape {
public:
	Shape(ReflectionType type, const Vec &color, const Vec &emission, bool isL = false) :
		reflType(type), c(color), e(emission), isLight(isL) {}
	virtual real Intersect(const Ray &r, Intersection *isect) const = 0;
	virtual bool Intersect(const Ray &r, Intersection *isect, real *t) const = 0;
	virtual bool Intersect(const Ray &r) const = 0;
	virtual Vec Sample(real *pdf, Vec rand) const = 0;
	virtual Vec Sample(const Intersection &isect, real *pdf, Vec u) const = 0;

	virtual real Pdf(const Intersection &isect, const Vec &wi) const {
		Ray ray(isect.hit, wi);
		Intersection lightPoint;
		real t;
		if (!Intersect(ray, &lightPoint, &t)) {
			return 0;
		}
		Vec d = lightPoint.hit - isect.hit;
		real pdf = d.dot(d) / (std::abs(lightPoint.n.dot(-1 * wi)) * Area());
		if (std::isinf(pdf)) return 0;
		return pdf;
	}

	virtual std::shared_ptr<BSDF> GetBSDF(const Intersection &isect, TransportMode mode = TransportMode::Radiance) const {
		if (reflType == DIFF) {
			return std::dynamic_pointer_cast<BSDF>(std::make_shared<DiffuseBSDF>(isect, c));
		}
		else if (reflType == SPEC) {
			return std::dynamic_pointer_cast<BSDF>(std::make_shared<SpecularBSDF>(isect, c));
		}
		else if (reflType == REFR) {
			return std::dynamic_pointer_cast<BSDF>(std::make_shared<TransmissionBSDF>(isect, c, mode));
		}
		return std::shared_ptr<BSDF>();
	}

	bool IsLight() const { return isLight; }

	virtual Vec GetNorm(const Vec &point) const = 0;

	virtual Vec GetEmission() const { return e; }

	int GetId() const { return shapeId; }

	virtual real Area() const = 0;

	friend class Scene;
private:
	ReflectionType reflType;
	Vec c, e;
	bool isLight;
	int shapeId;
};
