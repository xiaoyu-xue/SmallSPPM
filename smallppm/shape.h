#pragma once

#include "linagl.h"
#include "material.h"
#include "ray.h"
#include "intersection.h"
#include "bsdf.h"

NAMESPACE_BEGIN

class Shape {
public:
	virtual bool Intersect(const Ray &r, Intersection *isect, real *t) const = 0;
	virtual bool Intersect(const Ray &r) const = 0;
	virtual Intersection Sample(real *pdf, const Vec2 &rand) const = 0;
	virtual Intersection Sample(const Intersection &isect, real *pdf, const Vec2 &u) const = 0;

	virtual real Pdf(const Intersection &isect, const Vec3 &wi) const {
		//Ray ray(isect.hit, wi);
		Ray ray = isect.SpawnRay(wi);
		Intersection lightPoint;
		real t = Inf;
		//if (debugPixel == 1) {
		//	std::cout << "Pdf ray " << ray << std::endl;
		//}
		if (!Intersect(ray, &lightPoint, &t)) {
			//if (debugPixel == 1) {
			//	std::cout << "FindLightPdf_t: "<< t << std::endl;
			//}
			return 0;
		}
		Vec3 d = lightPoint.hit - isect.hit;
		real pdf = d.Dot(d) / (std::abs(lightPoint.n.Dot(-1 * wi)) * Area());
		//if (debugPixel == 1) {
		//	std::cout << "FindLightPdf_ pdf: " << pdf << std::endl;
		//}
		if (std::isinf(pdf)) return 0;
		return pdf;
	}

	virtual Vec3 GetNorm(const Vec3 &point) const = 0;

	int GetId() const { return shapeId; }

	virtual real Area() const = 0;

	friend class Scene;
private:
	int shapeId;
};

NAMESPACE_END