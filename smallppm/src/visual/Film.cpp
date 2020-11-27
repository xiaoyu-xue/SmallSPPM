#include "Film.h"
#include "image/ImageIO.h"

NAMESPACE_BEGIN

void Film::SaveImage() {
	WriteToPixelBuffer();
	ImageIO::WriteImage(mFilename, mImageBuffer, resX, resY);
}

NAMESPACE_END