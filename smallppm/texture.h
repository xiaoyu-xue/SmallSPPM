#pragma once
#include "def.h"
#include "linagl.h"
#include "intersection.h"

NAMESPACE_BEGIN

class Texture {
public:
	virtual ~Texture() { }

	virtual Vec3 Sample(const Vec3& coord) const {
		return Vec3(0, 0, 0);
	}

	virtual Vec3 Sample(const Vec2& coord) const {
		return Vec3(0, 0, 0);
	}

	virtual Vec3 Sample(const Intersection& isect) const {
		return Vec3(0, 0, 0);
	}
};

NAMESPACE_END