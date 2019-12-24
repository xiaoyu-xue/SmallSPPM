#pragma once

#include "camera.h"

NAMESPACE_BEGIN

class PinHoleCamera : public ProjectiveCamera {
public:
	PinHoleCamera(const std::shared_ptr<Film> &pFilm, const Vec3 &position,
		const Vec3 &pCz, const Vec3 &pCx, const Vec3 &pCy, real pFovy, real dis) :
		ProjectiveCamera(pFilm, position, pCz, pCx, pCy, pFovy, dis) {
		//Initialize();
	}

	Ray GenerateRay(int pixelX, int pixelY, const Vec2 &u, real offset) const override {
		//Vec3 dir = cx * ((pixelX + u.x) / film->resX - 0.5f) * film->width +
		//	cy * (-(pixelY + u.y) / film->resY + 0.5f) * film->heigh + cz * filmDistance;

		Vec3 sampledPos = RasterToWorld(Vec3(pixelX + u.x, pixelY + u.y, 0));

		Vec3 dir = (sampledPos - pos).Norm();
		return Ray(pos + dir * offset, dir, Inf, 0);
	}

	Vec3 We(const Ray &ray) const override {
		real pdfA = 1.0; // for the pinhole camera
		real area = film->area;
		real cosTheta = cz.Dot(ray.d);
		real cos2Theta = cosTheta * cosTheta;
		real value = filmDistance * filmDistance * pdfA / (area * cos2Theta * cos2Theta);
		return Vec3(value, value, value);
	}

	Vec3 Sample_Wi(const Intersection &isect, real *pdfW, Vec3 *wi, Vec3 u) const override {
		*wi = (pos - isect.hit);
		real distance = wi->Length();
		*wi = wi->Norm();
		real cosTheta = cz.Dot(-1 * (*wi));
		*pdfW = 1.f * (distance * distance) / cosTheta;
		//*PdfW = 1.0 * (dis / CosTheta) * (dis / CosTheta) / CosTheta;
		return We(Ray(isect.hit, -1 * (*wi)));
	}

	real PdfPos() const override {
		return 1.0;
	}

	real PdfDir(const Ray &cameraRay) const override {
		real filmArea = film->Area();
		real cosTheta = std::abs(cz.Dot(cameraRay.d));
		real cos2Theta = cosTheta * cosTheta;
		return filmDistance * filmDistance / (filmArea * cos2Theta * cosTheta);
	}
private:
	//Vec3 pos, cx, cy, cz;
	//real fovy;
	//real filmDistance;

};

NAMESPACE_END