#include "Film.h"
#include "image/ImageIO.h"

NAMESPACE_BEGIN

void Film::SaveImage() {
	WriteToPixelBuffer();
	ImageIO::WriteImage(filename, imageBuffer, resX, resY);
}

NAMESPACE_END