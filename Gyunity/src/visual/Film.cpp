#include "Film.h"
#include "image/ImageIO.h"

GY_NAMESPACE_BEGIN

void Film::SaveImage() {
	WriteToPixelBuffer();
	ImageIO::WriteImage(mFilename, mImageBuffer, resX, resY);
}

GY_NAMESPACE_END