#pragma once

#include "visual/Intersection.h"


GYT_NAMESPACE_BEGIN

class PathVertex
{
public:

	Intersection mIsect;
	Vec3 mThroughput;
	real mPdfFwd = 0.0, mPdfPrev = 0.0;
	bool mIsDelta = false;
};


GYT_NAMESPACE_END