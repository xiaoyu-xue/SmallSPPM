#include "intersection.h"
#include "primitive.h"

void Intersection::ComputeScatteringFunction(TransportMode mode) {
	primitive->ComputeScatteringFunction(this, mode);
}