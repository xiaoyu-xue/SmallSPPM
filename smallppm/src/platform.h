#pragma once

#ifdef _MSC_VER
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
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
#endif