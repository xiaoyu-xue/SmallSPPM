#pragma once
#include "common/Core.h"
#include "math/Ray.h"
#include "visual/Intersection.h"
#include "visual/Primitive.h"
#include "math/AABB.h"

GYT_NAMESPACE_BEGIN

class Accelerator {
public:
	virtual ~Accelerator(){}
	virtual bool Intersect(const Ray& r, Intersection* isect) const = 0;
	virtual bool Intersect(const Ray& r) const = 0;
	virtual void SetPrimitives(const std::vector<std::shared_ptr<Primitive>> &pPrimitives){}
	virtual AABB WorldBound() const = 0;
};

GYT_NAMESPACE_END
