#include "Intersection.h"
#include "Primitive.h"
#include "system/Memory.h"

NAMESPACE_BEGIN

void Intersection::ComputeScatteringFunction(MemoryPool &arena, TransportMode mode) {
	primitive->ComputeScatteringFunction(this, arena, mode);
}

NAMESPACE_END