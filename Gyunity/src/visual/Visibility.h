#pragma once

#include "Intersection.h"

GYT_NAMESPACE_BEGIN

class Scene;
class StateSequence;

class VisibilityTester {
private:
	Intersection mP0, mP1;
public:
	VisibilityTester() {}
	VisibilityTester(const Intersection &p0, const Intersection &p1)
		: mP0(p0), mP1(p1) 
	{
	}

	GYT_FORCE_INLINE const Intersection &P0() const { return mP0; }

	GYT_FORCE_INLINE const Intersection &P1() const { return mP1; }

	bool Unoccluded(const Scene& scene) const;

	Vec3 Tr(const Scene& scene, StateSequence& rand) const;

};

GYT_NAMESPACE_END