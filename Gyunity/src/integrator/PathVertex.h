#pragma once

#include "visual/Intersection.h"


GYT_NAMESPACE_BEGIN

class PathVertex
{
public:

	Intersection isect;
	Vec3 Throughput;
	double PdfFwd = 0.0, PdfPrev = 0.0;
	Camera* camera;
};


GYT_NAMESPACE_END