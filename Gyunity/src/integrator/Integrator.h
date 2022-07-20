#pragma once

#include "common/Core.h"
#include "core/Scene.h"
#include "core/Intersection.h"
#include "core/Visibility.h"
#include "common/DebugUtils.h"

GYT_NAMESPACE_BEGIN

GYT_FORCE_INLINE real PowerHeuristic(int nf, real fPdf, int ng, real gPdf) 
{
	real f = nf * fPdf, g = ng * gPdf;
	return (f * f) / (f * f + g * g);
}

class Integrator
{
public:
	virtual void Render(const Scene &scene, const Camera &camera) = 0;

	virtual ~Integrator() {}

	static Vec3 DirectIllumination(const Scene& scene, const Intersection& isect, real uLight, 
		const Vec2& u, const Vec3& v, StateSequence& rand, bool handleMedia = false);

};

GYT_NAMESPACE_END