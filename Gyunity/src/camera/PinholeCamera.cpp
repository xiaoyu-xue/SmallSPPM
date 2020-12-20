#include "PinholeCamera.h"
#include "visual/Film.h"

GYT_NAMESPACE_BEGIN

Vec3 PinHoleCamera::We(const Ray& ray) const  {
	real pdfA = 1.0; // for the pinhole camera
	real area = film->mArea;
	real cosTheta = cz.Dot(ray.mDir);
	real cos2Theta = cosTheta * cosTheta;
	real value = filmDistance * filmDistance * pdfA / (area * cos2Theta * cos2Theta);
	return Vec3(value, value, value);
}

real PinHoleCamera::PdfDir(const Ray& cameraRay) const {
	real filmArea = film->Area();
	real cosTheta = std::abs(cz.Dot(cameraRay.mDir));
	real cos2Theta = cosTheta * cosTheta;
	return filmDistance * filmDistance / (filmArea * cos2Theta * cosTheta);
}

GYT_NAMESPACE_END