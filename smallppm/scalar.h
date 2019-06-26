#pragma once

NAMESPACE_BEGIN

template <class T1, class T2, class T3>
T1 Clamp(const T1& tVal, const T2& tMin, const T3& max)
{
	if (tVal < tMin) return tMin;
	if (tVal > max) return max;
	return tVal;
}

NAMESPACE_END