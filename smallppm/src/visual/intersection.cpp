#include "intersection.h"
#include "primitive.h"
#include "system/memory.h"

void Intersection::ComputeScatteringFunction(MemoryArena &arena, TransportMode mode) {
	primitive->ComputeScatteringFunction(this, arena, mode);
}