#pragma once
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#define USING_DOUBLE


#include <stdint.h>
#include <algorithm>
#include "def.h"

NAMESPACE_BEGIN

/***************** Type ******************/
using uchar = unsigned char;

using int8 = int8_t;
using uint8 = uint8_t;

using int16 = int16_t;
using uint16 = uint16_t;

using int32 = int32_t;
using uint32 = uint32_t;

using int64 = int64_t;
using uint64 = uint64_t;

using float32 = float;
using float64 = double;

#ifdef USING_DOUBLE
using real = double;
using real_bit = uint64;
#else
using real = float;
using real_bit = uint32;
#endif

#ifdef _WIN64
#define FORCE_INLINE __forceinline
#else
#define FORCE_INLINE inline __attribute__((always_inline))
#endif




/****************** Constant *********************/

constexpr real PI = 3.14159265358979;
constexpr real INV_PI = 0.31830988618379067154;
constexpr real PiOver2 = 1.57079632679489661923;
constexpr real PiOver4 = 0.78539816339744830961;
constexpr real eps = 1e-6;
constexpr real Inf = 1e20;
constexpr real rayeps = 1e-4;
constexpr real shadowRayEps = 1e-4;
constexpr real nEps = 1e-6;
constexpr float64 NumericalEps = 1e-6;
constexpr real MachineEps = std::numeric_limits<real>::epsilon() * 0.5;
constexpr real MaxReal = std::numeric_limits<real>::max();
constexpr real Infinity = std::numeric_limits<real>::infinity();

NAMESPACE_END