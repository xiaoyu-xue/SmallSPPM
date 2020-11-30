#pragma once


#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

// Platforms

// Windows
#if defined(_WIN64)
#define GYT_PLATFORM_WINDOWS
#endif

#if defined(_WIN32) && !defined(_WIN64)
static_assert(false, "32-bit Windows systems are not supported")
#endif

// Compilers
// MSVC
#if defined(_MSC_VER)
#define GYT_COMPILER_MSVC
#endif

#   pragma warning(disable : 4275)
#   pragma warning(disable : 4267)
#   pragma warning(disable : 4251) // 'field' : class 'A' needs to have dll-interface to be used by clients of class 'B'
#   pragma warning(disable : 4800) // 'type' : forcing value to bool 'true' or 'false' (performance warning)
#   pragma warning(disable : 4996) // Secure SCL warnings
#   pragma warning(disable : 5030)
#   pragma warning(disable : 4324)
#   pragma warning(disable : 4201)
#   pragma warning(disable : 4305)
#   pragma warning(disable : 4244)
#   pragma warning(disable : 4805)
#   pragma warning(disable : 4018)