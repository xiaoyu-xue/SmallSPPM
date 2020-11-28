#pragma once
#include "common/Core.h"
#include "math/Linagl.h"
#include "visual/Intersection.h"

GY_NAMESPACE_BEGIN

template<typename T>
class Texture {
public:
	Texture() = default;

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

GY_NAMESPACE_END