#include "film.h"
#include "image/image_io.h"

NAMESPACE_BEGIN

void Film::SaveImage() {
	WriteToPixelBuffer();
	ImageIO::WriteImage(filename, imageBuffer, resX, resY);
}

NAMESPACE_END