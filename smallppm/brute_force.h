#pragma once
#include "def.h"
#include "accelerator.h"
#include "primitive.h"

NAMESPACE_BEGIN

class BruteForce : public Accelerator {
public:

	BruteForce() { }

	BruteForce(const std::vector<std::shared_ptr<Primitive>> primitives) :
		primitives(primitives) {

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
private:
	std::vector<std::shared_ptr<Primitive>> primitives;
};

NAMESPACE_END