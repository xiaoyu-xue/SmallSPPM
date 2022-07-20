#pragma once

#include "core/Intersection.h"
#include "ForwardDecl.h"

GYT_NAMESPACE_BEGIN

class PathVertex
{
public:

	Intersection mIsect;
	Vec3 mThroughput;
	real mPdfFwd = 0.0, mPdfPrev = 0.0;
	bool mIsDelta = false;
	Camera* mpCamera;
	Light* mpLight;
};


GYT_NAMESPACE_END