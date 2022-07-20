#pragma once
#include "common/Core.h"
#include "accelerator/Accelerator.h"
#include "core/Primitive.h"

GYT_NAMESPACE_BEGIN

class BruteForce : public Accelerator {
public:

	BruteForce() { }

	BruteForce(const std::vector<std::shared_ptr<Primitive>> primitives) :
		primitives(primitives) {
		for (int i = 0; i < primitives.size(); ++i) {
			bounds = Union(bounds, primitives[i]->WorldBound());
		}
	}

	bool Intersect(const Ray& r, Intersection* isect) const override {
		int n = (int)primitives.size();
		bool intersected = false;
		for (int i = 0; i < n; ++i) {
			if (primitives[i]->Intersect(r, isect)) {
				intersected = true;
			}
		}
		return  intersected;
	}

	bool Intersect(const Ray& r) const override {
		for (auto primitive : primitives) {
			if (primitive->Intersect(r)) {
				return true;
			}
		}
		return false;
	}

	void SetPrimitives(const std::vector<std::shared_ptr<Primitive>> &pPrimitives) override {
		primitives = pPrimitives;
	}

	AABB WorldBound() const { return bounds; }
private:
	std::vector<std::shared_ptr<Primitive>> primitives;
	AABB bounds;
};

GYT_NAMESPACE_END