#pragma once
#include "def.h"
#include "math/ray.h"
#include "visual/intersection.h"
#include "visual/primitive.h"
#include "math/AABB.h"

NAMESPACE_BEGIN

class Accelerator {
public:
	virtual ~Accelerator(){}
	virtual bool Intersect(const Ray& r, Intersection* isect) const = 0;
	virtual bool Intersect(const Ray& r) const = 0;
	virtual void SetPrimitives(const std::vector<std::shared_ptr<Primitive>> &pPrimitives){}
	virtual AABB WorldBound() const = 0;
};

NAMESPACE_END
