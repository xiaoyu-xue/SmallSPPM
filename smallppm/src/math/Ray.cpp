#include "Ray.h"

NAMESPACE_BEGIN

std::ostream& operator<<(std::ostream &os, const Ray &ray) {
	os << "ray: " << "orig: " << ray.mOrig << " dir: " << ray.mDir << " tMin: " << ray.m_tMin << " tMax: " << ray.m_tMax;
	return os;
}

NAMESPACE_END