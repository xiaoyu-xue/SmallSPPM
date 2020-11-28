#pragma once

#include "math/Linagl.h"
#include "math/Ray.h"
#include "visual/Intersection.h"
#include "math/Transform.h"
#include "visual/Film.h"

GY_NAMESPACE_BEGIN


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
		const Vec3& pCz, const Vec3& pCx, const Vec3& pCy, real pFovy, real dis);
protected:
	void Initialize();
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
GY_NAMESPACE_END