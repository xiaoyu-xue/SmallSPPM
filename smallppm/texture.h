#pragma once
#include "def.h"
#include "linagl.h"
#include "intersection.h"

NAMESPACE_BEGIN

template<typename T>
class Texture {
public:
	virtual ~Texture() { }

	virtual T Sample(const Vec3& coord) const = 0;

	virtual T Sample(const Vec2& coord) const = 0;

	virtual T Sample(const Intersection& isect) const = 0;
};

NAMESPACE_END