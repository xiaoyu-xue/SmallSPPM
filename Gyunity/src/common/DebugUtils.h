#pragma once
#include "common/Core.h"


GYT_NAMESPACE_BEGIN

extern int debugPixel[100];

#define DEBUG_PIXEL(pixelX, pixelY, threadIndex)\
if (x == (pixelX) && y == pixelY){ \
	debugPixel[(threadIndex)] = 1;				   \
}								   \
else{                              \
	debugPixel[(threadIndex)] = 0;				   \
}								   \


#define DEBUG_PIXEL_IF(threadIndex)\
if(debugPixel[(threadIndex)] == 1)		\

GYT_NAMESPACE_END