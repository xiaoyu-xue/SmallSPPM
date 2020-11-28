#include "Intersection.h"
#include "Primitive.h"
#include "system/Memory.h"

GY_NAMESPACE_BEGIN

void Intersection::ComputeScatteringFunction(MemoryPool &arena, TransportMode mode) {
	primitive->ComputeScatteringFunction(this, arena, mode);
}

GY_NAMESPACE_END