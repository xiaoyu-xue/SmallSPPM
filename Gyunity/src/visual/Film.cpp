#include "Film.h"
#include "image/ImageIO.h"

GYT_NAMESPACE_BEGIN

void Film::SaveImage() {
	WriteToPixelBuffer();
	ImageIO::WriteImage(mFilename, mImageBuffer, resX, resY);
}

GYT_NAMESPACE_END