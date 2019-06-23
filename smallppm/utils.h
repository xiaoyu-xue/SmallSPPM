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