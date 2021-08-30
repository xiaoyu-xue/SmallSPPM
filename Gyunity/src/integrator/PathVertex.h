#pragma once

#include "visual/Intersection.h"


GYT_NAMESPACE_BEGIN

class PathVertex
{
public:

	Intersection mIsect;
	Vec3 mThroughput;
	double mPdfFwd = 0.0, mPdfPrev = 0.0;
	Camera* mpCamera;
};


GYT_NAMESPACE_END