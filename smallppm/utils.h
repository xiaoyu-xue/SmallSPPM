#pragma once
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif
#define USING_DOUBLE


#include <stdint.h>


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
using real = float64;
#else
using real = float32;
#endif

#ifdef _WIN64
#define FORCE_INLINE __forceinline
#else
#define FORCE_INLINE inline __attribute__((always_inline))
#endif




/****************** Constant *********************/

const real PI = 3.14159265358979;
const real INV_PI = 0.31830988618379067154;
const real PiOver2 = 1.57079632679489661923;
const real PiOver4 = 0.78539816339744830961;
const real eps = 1e-6;
const real Inf = 1e20;
const real rayeps = 1e-4;
const real shadowRayEps = 1e-4;
const real nEps = 1e-6;