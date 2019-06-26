#pragma once
#include "def.h"
#include "utils.h"

NAMESPACE_BEGIN

inline uint32 reverse_bit32(uint32 n) {
	n = (n << 16) | (n >> 16);
	n = ((n & 0x00ff00ff) << 8) | ((n & 0xff00ff00) >> 8);
	n = ((n & 0x0f0f0f0f) << 4) | ((n & 0xf0f0f0f0) >> 4);
	n = ((n & 0x33333333) << 2) | ((n & 0xcccccccc) >> 2);
	n = ((n & 0x55555555) << 1) | ((n & 0xaaaaaaaa) >> 1);
	return n;
}

inline uint64 reverse_bit64(uint64 n) {
	uint64 n0 = reverse_bit32((uint32)n);
	uint64 n1 = reverse_bit32((uint32)(n >> 32));
	return (n0 << 32) | n1;
}

int IsPrime(int a) noexcept;

NAMESPACE_END