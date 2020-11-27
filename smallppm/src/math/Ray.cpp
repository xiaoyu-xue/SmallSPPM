#include "Ray.h"

NAMESPACE_BEGIN

std::ostream& operator<<(std::ostream &os, const Ray &ray) {
	os << "ray: " << "orig: " << ray.mOrig << " dir: " << ray.mDir << " tMin: " << ray.tMin << " tMax: " << ray.tMax;
	return os;
}

NAMESPACE_END