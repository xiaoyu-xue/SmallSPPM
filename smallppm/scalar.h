#pragma once
#include <algorithm>
#include <utility>

NAMESPACE_BEGIN

template <class T1, class T2, class T3>
T1 Clamp(const T1& tVal, const T2& tMin, const T3& max)
{
	if (tVal < tMin) return tMin;
	if (tVal > max) return max;
	return tVal;
}

FORCE_INLINE real Radians(real deg) { return (PI / 180) * deg; }

FORCE_INLINE real Degrees(real rad) { return (180 / PI) * rad; }

NAMESPACE_END