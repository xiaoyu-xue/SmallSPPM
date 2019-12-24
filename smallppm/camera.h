#pragma once

#include "film.h"
#include "linagl.h"
#include "ray.h"
#include "intersection.h"
#include "transform.h"

NAMESPACE_BEGIN

class Camera {
public:
	Camera(const std::shared_ptr<Film> &pFilm) {
		film = pFilm;
	}
	virtual Ray GenerateRay(int pixelX, int pixelY, const Vec2 &u, real offset = 0.0) const = 0;
	virtual Vec3 We(const Ray &ray) const = 0;
	virtual Vec3 Sample_Wi(const Intersection &isect, real *pdfW, Vec3 *wi, Vec3 u) const = 0;
	virtual real PdfPos() const = 0;
	virtual real PdfDir(const Ray &cameraRay) const = 0;
	virtual std::shared_ptr<Film> GetFilm() const { return film; }
protected:

	std::shared_ptr<Film> film;
};

class ProjectiveCamera : public Camera {
public:
	ProjectiveCamera(const std::shared_ptr<Film>& pFilm, const Vec3& position,
		const Vec3& pCz, const Vec3& pCx, const Vec3& pCy, real pFovy, real dis) :
		Camera(pFilm), pos(position), cx(pCx), cy(pCy), cz(pCz), fovy(pFovy), filmDistance(dis) {

		Initialize();

		CameraToWorld = 
			Transform(cx.x, cy.x, cz.x, pos.x,
					  cx.y, cy.y, cz.y, pos.y,
					  cx.z, cy.z, cz.z, pos.z,
					  0, 0, 0, 1.f);

		WorldToCamera = Inverse(CameraToWorld);

		CameraToNDC = Transform::Perspective(fovy, film->aspect, filmDistance, filmDistance, 100000.f);

		NDCToCamera = Inverse(CameraToNDC);

		NDCToRaster =
			Transform::Scale(real(film->resX), real(film->resY), 1) *
			Transform::Scale(0.5f, -0.5f, 1) * 
			Transform::Translate(Vec3(1, -1, 0));

		RasterToNDC = Inverse(NDCToRaster);

		CameraToRaster = NDCToRaster * CameraToNDC;

		RasterToCamera = NDCToCamera * RasterToNDC;

		RasterToWorld = CameraToWorld * RasterToCamera;

	}
protected:

	void Initialize() {
		Vec3 filmCenter = pos + czz * filmDistance;
		//std::cout << "film cetner: " << filmCenter << std::endl;
		film->heigh = filmDistance * std::tan(fovy * 0.5f * PI / 180) * 2.f;
		film->width = film->heigh * film->aspect;
		film->area = film->width * film->heigh;

		film->LU = filmCenter + cy * film->heigh * 0.5 - cx * film->width * 0.5;
		film->LL = filmCenter - cy * film->heigh * 0.5 - cx * film->width * 0.5;
		film->RU = filmCenter + cy * film->heigh * 0.5 + cx * film->width * 0.5;
		film->RL = filmCenter - cy * film->heigh * 0.5 + cx * film->width * 0.5;
		std::cout << "LL: " << film->LL << std::endl
			<< "LU: " << film->LU << std::endl
			<< "RL: " << film->RL << std::endl
			<< "RR: " << film->RU << std::endl;
	}
public:
	Vec3 pos, cx, cy, cz, czz;
	real fovy;
	real filmDistance;

	Transform CameraToWorld, WorldToCamera;
	Transform CameraToNDC, NDCToCamera;
	Transform NDCToRaster, RasterToNDC;
	Transform CameraToRaster, RasterToCamera;
	Transform RasterToWorld;
};
NAMESPACE_END