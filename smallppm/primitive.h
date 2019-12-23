#pragma once
#include "intersection.h"
class Primitive {
public:

	virtual bool Intersect(const Ray& r, Intersection* isect, real* t) const = 0;
	virtual bool Intersect(const Ray& r) const = 0;

	virtual ~Primitive(){}
};