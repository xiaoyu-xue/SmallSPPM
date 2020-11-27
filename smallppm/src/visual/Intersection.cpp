#include "Intersection.h"
#include "Primitive.h"
#include "system/Memory.h"

void Intersection::ComputeScatteringFunction(MemoryPool &arena, TransportMode mode) {
	primitive->ComputeScatteringFunction(this, arena, mode);
}