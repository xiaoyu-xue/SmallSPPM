#include "Ray.h"

NAMESPACE_BEGIN

std::ostream& operator<<(std::ostream &os, const Ray &ray) {
	os << "ray: " << "orig: " << ray.o << " dir: " << ray.d << " tMin: " << ray.tMin << " tMax: " << ray.tMax;
	return os;
}

NAMESPACE_END