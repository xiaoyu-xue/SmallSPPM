#include "Intersection.h"
#include "Primitive.h"
#include "system/Memory.h"

GYT_NAMESPACE_BEGIN

void Intersection::ComputeScatteringFunction(MemoryPool &arena, TransportMode mode) {
	primitive->ComputeScatteringFunction(this, arena, mode);
}

GYT_NAMESPACE_END