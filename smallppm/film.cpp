#include "film.h"
#include "image_io.h"

NAMESPACE_BEGIN

void Film::SaveImage() {
	ImageIO::WriteImage(filename, imageBuffer, resX, resY);
}

NAMESPACE_END