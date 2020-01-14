#pragma once
#include "utils.h"


NAMESPACE_BEGIN

extern int debugPixel;

#define DEBUG_PIXEL(pixelX, pixelY)\
if (x == (pixelX) && y == pixelY){ \
	debugPixel = 1;				   \
}								   \
else{                              \
	debugPixel = 0;				   \
}								   \


#define DEBUG_PIXEL_IF()\
if(debugPixel == 1)		\

NAMESPACE_END