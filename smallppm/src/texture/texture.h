#pragma once
#include "def.h"
#include "math/linagl.h"
#include "visual/intersection.h"

NAMESPACE_BEGIN

template<typename T>
class Texture {
public:
	virtual ~Texture() { }

	virtual T Sample(const Vec3& coord) const {
		return T(0);
	}

	virtual T Sample(const Vec2& coord) const {
		return T(0);
	}

	virtual T Sample(const Intersection& isect) const {
		return T(0);
	}
};

NAMESPACE_END