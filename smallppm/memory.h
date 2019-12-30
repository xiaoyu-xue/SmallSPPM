#pragma once
#include "def.h"

#define L1_CACHE_LINE_SIZE 64

NAMESPACE_BEGIN

void* AllocAligned(size_t size);
template <typename T>
T* AllocAligned(size_t count) {
	return (T*)AllocAligned(count * sizeof(T));
}

void FreeAligned(void*);

class MemoryArena {

};



NAMESPACE_END