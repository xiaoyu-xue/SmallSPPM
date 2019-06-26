#pragma once

#include "camera.h"

NAMESPACE_BEGIN

class PinHoleCamera : public Camera {
public:
	PinHoleCamera(const std::shared_ptr<Film> &pFilm, const Vec &position,
		const Vec &pCz, const Vec &pCx, const Vec &pCy, real pFovy, real dis) :
		Camera(pFilm), pos(position), cz(pCz), cx(pCx), cy(pCy), fovy(pFovy), filmDistance(dis) {
		Initialize();
	}

	void Initialize() {
		Vec filmCenter = pos + cz * filmDistance;
		//std::cout << "film cetner: " << filmCenter << std::endl;
		film->heigh = filmDistance * std::tan(fovy * 0.5 * PI / 180) * 2.0;
		film->width = film->heigh * film->aspect;
		film->area = film->width * film->heigh;
		/*
		film->LU = filmCenter + cy * film->heigh * 0.5 - cx * film->width * 0.5;
		film->LL = filmCenter - cy * film->heigh * 0.5 - cx * film->width * 0.5;
		film->RU = filmCenter + cy * film->heigh * 0.5 + cx * film->width * 0.5;
		film->RL = filmCenter - cy * film->heigh * 0.5 + cx * film->width * 0.5;
		*/
	}

	Ray GenerateRay(int pixelX, int pixelY, Vec u, real offset) const {

		Vec dir = cx * ((pixelX + u.x) / film->resX - 0.5) * film->width +
			cy * (-(pixelY + u.y) / film->resY + 0.5) * film->heigh + cz * filmDistance;
		dir = dir.norm();

		return Ray(pos + dir * offset, dir);
	}

	Vec We(const Ray &ray) const {
		real pdfA = 1.0; // for the pinhole camera
		real area = film->area;
		real cosTheta = cz.dot(ray.d);
		real cos2Theta = cosTheta * cosTheta;
		real value = filmDistance * filmDistance * pdfA / (area * cos2Theta * cos2Theta);
		return Vec(value, value, value);
	}

	Vec Sample_Wi(const Intersection &isect, real *pdfW, Vec *wi, Vec u) const {
		*wi = (pos - isect.hit);
		real distance = wi->length();
		*wi = wi->norm();
		real cosTheta = cz.dot(-1 * (*wi));
		*pdfW = 1.0 * (distance * distance) / cosTheta;
		//*PdfW = 1.0 * (dis / CosTheta) * (dis / CosTheta) / CosTheta;
		return We(Ray(isect.hit, -1 * (*wi)));
	}

	real PdfPos() const {
		return 1.0;
	}

	real PdfDir(const Ray &cameraRay) const {
		real filmArea = film->Area();
		real cosTheta = std::abs(cz.dot(cameraRay.d));
		real cos2Theta = cosTheta * cosTheta;
		return filmDistance * filmDistance / (filmArea * cos2Theta * cosTheta);
	}
private:
	Vec pos, cx, cy, cz;
	real fovy;
	real filmDistance;

};

NAMESPACE_END