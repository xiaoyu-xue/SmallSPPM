#pragma once

#include "linagl.h"
#include "material.h"
#include "ray.h"
#include "intersection.h"
#include "bsdf.h"
#include "transform.h"

NAMESPACE_BEGIN

class Shape {
public:
	Shape() {}
	Shape(const Transform* ObjectToWorld, const Transform* WorldToObject) :
		ObjectToWorld(ObjectToWorld), WorldToObject(WorldToObject) {

	}
	virtual ~Shape(){ }
	virtual bool Intersect(const Ray &r, Intersection *isect, real *t) const = 0;
	virtual bool Intersect(const Ray &r) const = 0;
	virtual Intersection Sample(real *pdf, const Vec2 &rand) const = 0;
	virtual Intersection Sample(const Intersection& isect, real* pdf, const Vec2& u) const;
	virtual AABB ObjectBound() const = 0;
	virtual AABB WorldBould() const{
		if(ObjectToWorld) return (*ObjectToWorld)(ObjectBound());
		return ObjectBound();
	}
	virtual real Pdf(const Intersection &isect, const Vec3 &wi) const {
		//Ray ray(isect.hit, wi);
		Ray ray = isect.SpawnRay(wi);
		Intersection lightPoint;
		real t = Inf;
		if (!Intersect(ray, &lightPoint, &t)) {
			return 0;
		}
		Vec3 d = lightPoint.hit - isect.hit;
		real pdf = d.Dot(d) / (std::abs(lightPoint.n.Dot(-1 * wi)) * Area());
		if (std::isinf(pdf)) return 0;
		return pdf;
	}

	virtual Vec3 GetNorm(const Vec3 &point) const = 0;

	int64 GetId() const { return shapeId; }

	virtual real Area() const = 0;

	friend class Scene;

	int64 shapeId;
protected:

	const Transform* ObjectToWorld, * WorldToObject;
};

NAMESPACE_END