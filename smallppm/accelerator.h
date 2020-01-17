#pragma once
#include "def.h"
#include "ray.h"
#include "intersection.h"
#include "primitive.h"

NAMESPACE_BEGIN

class Accelerator {
public:
	virtual ~Accelerator(){}
	virtual bool Intersect(const Ray& r, Intersection* isect) const = 0;
	virtual bool Intersect(const Ray& r) const = 0;
	virtual void SetPrimitives(const std::vector<std::shared_ptr<Primitive>> &pPrimitives){}
};

NAMESPACE_END
