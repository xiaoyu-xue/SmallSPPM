#pragma once

#include "intersection.h"

NAMESPACE_BEGIN

class VisibilityTester {
public:
	VisibilityTester() {}
	VisibilityTester(const Intersection &p0, const Intersection &p1)
		: p0(p0), p1(p1) {}
	const Intersection &P0() const { return p0; }
	const Intersection &P1() const { return p1; }
	bool Unoccluded(const Scene &scene) const {
		Ray ray = p0.SpawnTo(p1);
		Intersection isect;
		real t;
		//DEBUG_PIXEL_IF() {
		//	std::cout << "p0: " << p0.shapeId << " " << "p1: " << p1.shapeId << std::endl;
		//	std::cout << scene.Intersect(ray, &isect, &t) << " " << isect.shapeId << std::endl;
		//}
		//return scene.Intersect(ray, &isect, &t) && isect.shapeId != p1.shapeId;
		return scene.Intersect(ray);
	};

private:
	Intersection p0, p1;
};

NAMESPACE_END