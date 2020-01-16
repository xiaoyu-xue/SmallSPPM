#include "intersection.h"
#include "primitive.h"

void Intersection::ComputeScatteringFunction(MemoryArena &arena, TransportMode mode) {
	primitive->ComputeScatteringFunction(this, arena, mode);
}