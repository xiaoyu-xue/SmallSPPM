#pragma once

#include "common/Core.h"
#include "math/Linagl.h"
#include "math/Ray.h"
#include "visual/Intersection.h"
#include "math/Transform.h"
#include "visual/Film.h"

GYT_NAMESPACE_BEGIN


class Camera 
{
protected:
	std::shared_ptr<Film> mpFilm;
public:
	Vec3 mPos, mLookAt, mUp;
	Vec3 mCx, mCy, mCz;
	real mFovy;
	real mFilmDistance;

	Transform CameraToWorld, WorldToCamera;
	Transform CameraToNDC, NDCToCamera;
	Transform NDCToRaster, RasterToNDC;
	Transform CameraToRaster, RasterToCamera;
	Transform RasterToWorld;
public:
	Camera(std::shared_ptr<Film> pFilm) 
	{
		mpFilm = pFilm;
	}

	Camera(std::shared_ptr<Film> pFilm, const Vec3& position, const Vec3& cz, const Vec3& cx, const Vec3& cy, real fovy, real filmDis) :
		mPos(position), mCx(cx), mCy(cy), mCz(cz), mFovy(fovy), mFilmDistance(filmDis) {
		mpFilm = pFilm;
	}
	Camera(std::shared_ptr<Film> pFilm, const Vec3& position, const Vec3& lookAt, const Vec3& up, real fovy, real filmDis) :
		mPos(position), mLookAt(lookAt), mUp(up), mFovy(fovy), mFilmDistance(filmDis) {
		mpFilm = pFilm;
	}
	virtual Ray GenerateRay(int pixelX, int pixelY, const Vec2 &u, real offset = 0.0) const = 0;
	virtual Vec3 We(const Ray &ray) const = 0;
	virtual Vec3 Sample_Wi(const Intersection &isect, real *pdfW, Vec3 *wi, Vec3 u) const = 0;
	virtual real PdfPos() const = 0;
	virtual real PdfDir(const Ray &cameraRay) const = 0;
	virtual std::shared_ptr<Film> GetFilm() const { return mpFilm; }
};

class ProjectiveCamera : public Camera 
{
public:
	//Vec3 mPos, mLookAt, mUp;
	//Vec3 mCx, mCy, mCz;
	//real mFovy;
	//real mFilmDistance;

	//Transform CameraToWorld, WorldToCamera;
	//Transform CameraToNDC, NDCToCamera;
	//Transform NDCToRaster, RasterToNDC;
	//Transform CameraToRaster, RasterToCamera;
	//Transform RasterToWorld;
public:
	ProjectiveCamera(std::shared_ptr<Film> pFilm, const Vec3& position, const Vec3& cz, const Vec3& cx, const Vec3& cy, real fovy, real filmDis);
	ProjectiveCamera(std::shared_ptr<Film> pFilm, const Vec3& position, const Vec3& lookAt, const Vec3& up, real fovy, real filmDis);
	virtual ~ProjectiveCamera() { }
protected:
	void Initialize();
};

GYT_NAMESPACE_END