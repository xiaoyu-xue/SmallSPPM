#pragma once

#include "intersection.h"

NAMESPACE_BEGIN

class Scene;
class StateSequence;

class VisibilityTester {
public:
	VisibilityTester() {}
	VisibilityTester(const Intersection &p0, const Intersection &p1)
		: p0(p0), p1(p1) {}
	const Intersection &P0() const { return p0; }
	const Intersection &P1() const { return p1; }
	bool Unoccluded(const Scene& scene) const;
	Vec3 Tr(const Scene& scene, StateSequence& rand) const;

private:
	Intersection p0, p1;
};

NAMESPACE_END