#pragma once
#include "utils.h"
#include <algorithm>
#include "linagl.h"

NAMESPACE_BEGIN

template <typename T>
FORCE_INLINE
typename std::enable_if_t<std::is_floating_point<T>::value, bool>
Equal(const T &a, const T &b, float64 tolerance = NumericalEps) {
	return std::abs(a - b) <= tolerance;
}


template<typename T, typename TensorType TType = T::tensorType,
std::enable_if_t<TType == TensorType::VECTOR, int> = 0>
inline typename T::type Maximum(const T &t) {
	typename T::type ret = t(0);
	for (int i = 1; i < T::dim; ++i) {
		if (std::abs(ret) < std::abs(t(i))) {
			ret = t(i);
		}
	}
	return ret;
}

template<typename T, typename TensorType TType = T::tensorType,
	std::enable_if_t<TType == TensorType::MATRIX, int> = 0>
	inline typename T::type Maximum(const T &t) {
	typename T::type ret = t(0, 0);
	for (int i = 0; i < T::dim; ++i) {
		for (int j = 0; j < T::dim; ++j) {
			if (std::abs(ret) < std::abs(t(i, j))) {
				ret = t(i, j);
			}
		}
	}
	return ret;
}

template <typename T>
FORCE_INLINE
typename std::enable_if_t<!std::is_floating_point<T>::value, bool>
Equal(const T &a, const T &b, float64 tolerance = NumericalEps) {
	return std::abs(Maximum(a - b)) <= tolerance;
}

NAMESPACE_END