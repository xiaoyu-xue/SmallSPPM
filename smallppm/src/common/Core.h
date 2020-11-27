#pragma once

#include "Platform.h"
#include "Def.h"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <type_traits>
#include <cstdint>
#include <algorithm>
#include <memory>
#include <csignal>
#include <vector>
#include <list>

#define ISE_SSE

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

#if defined(PLATFORM_WINDOWS)
#define FORCE_INLINE __forceinline
#else
#define FORCE_INLINE inline __attribute__((always_inline))
#endif

#if defined(COMPILER_MSVC)
#define ALIGNED(x) __declspec(align(x))
#endif


/****************** Constant *********************/

constexpr real PI = (real)3.14159265358979;
constexpr real INV_PI = (real)0.31830988618379067154;
constexpr real INV_2PI = (real)0.15915494309189533577;
constexpr real INV_4PI = (real)0.07957747154594766788;
constexpr real PiOver2 = (real)1.57079632679489661923;
constexpr real PiOver4 = (real)0.78539816339744830961;
constexpr real Eps = (real)1e-6;
constexpr real Inf = std::numeric_limits<real>::infinity();
constexpr real RayEps = (real)1e-4;
constexpr real ShadowRayEps = (real)1e-4;
constexpr real PhtotonEdgeEps = (real)0.0009;
constexpr real NumericalEps = 1e-6;
constexpr real MachineEps = std::numeric_limits<real>::epsilon() * (real)0.5;
constexpr real MaxReal = std::numeric_limits<real>::max();
constexpr real Infinity = std::numeric_limits<real>::infinity();

NAMESPACE_END