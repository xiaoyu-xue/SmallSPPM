#include "linagl.h"

Vec operator*(real a, Vec b) { return Vec(a * b.x, a * b.y, a * b.z); }

std::ostream& operator<<(std::ostream &os, const Vec &v) {
	os << v.x << " " << v.y << " " << v.z;
	return os;
}