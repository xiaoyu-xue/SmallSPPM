#pragma once
#include <algorithm>
#include <utility>
#include <cmath>
#include "utils.h"

NAMESPACE_BEGIN

template <class T1, class T2, class T3>
T1 Clamp(const T1& tVal, const T2& tMin, const T3& max)
{
	if (tVal < tMin) return tMin;
	if (tVal > max) return max;
	return tVal;
}

template<typename T, typename V>
V Lerp(T a, V x0, V x1) {
	return (T(1) - a) * x0 + a * x1;
}

FORCE_INLINE real Radians(real deg) { return (PI / 180) * deg; }

FORCE_INLINE real Degrees(real rad) { return (180 / PI) * rad; }


FORCE_INLINE real Log2(real x) {
	const real invLog2 = 1.442695040888963387004650940071;
	return std::log(x) * invLog2;
}

FORCE_INLINE int Log2Int(uint32_t v) {
#if defined(_MSC_VER)
	unsigned long lz = 0;
	if (_BitScanReverse(&lz, v)) return lz;
	return 0;
#else
	return 31 - __builtin_clz(v);
#endif
}

FORCE_INLINE int Log2Int(int32_t v) { return Log2Int((uint32_t)v); }

FORCE_INLINE int Log2Int(uint64_t v) {
#if defined(_MSC_VER)
	unsigned long lz = 0;
#if defined(_WIN64)
	_BitScanReverse64(&lz, v);
#else
	if (_BitScanReverse(&lz, v >> 32))
		lz += 32;
	else
		_BitScanReverse(&lz, v & 0xffffffff);
#endif // _WIN64
	return lz;
#else  // IS_MSVC
	return 63 - __builtin_clzll(v);
#endif
}

FORCE_INLINE int Log2Int(int64_t v) { return Log2Int((uint64_t)v); }

template <typename T>
FORCE_INLINE  bool IsPowerOf2(T v) {
	return v && !(v & (v - 1));
}

inline int32_t RoundUpPow2(int32_t v) {
	v--;
	v |= v >> 1;
	v |= v >> 2;
	v |= v >> 4;
	v |= v >> 8;
	v |= v >> 16;
	return v + 1;
}

inline int64_t RoundUpPow2(int64_t v) {
	v--;
	v |= v >> 1;
	v |= v >> 2;
	v |= v >> 4;
	v |= v >> 8;
	v |= v >> 16;
	v |= v >> 32;
	return v + 1;
}

inline real Sgn(real val) {
	return (val > 0) ? 1.f : -1.f;
}


inline real SafeSqrt(real value) {
	return std::sqrt(std::max((real)0.0f, value));
}
NAMESPACE_END