#pragma once

#include "Camera.h"

GYT_NAMESPACE_BEGIN

class PinHoleCamera : public ProjectiveCamera 
{
public:
	PinHoleCamera(std::shared_ptr<Film> pFilm, const Vec3 &position, const Vec3 &cz, const Vec3 &cx, const Vec3 &cy, real fovy, real dis) 
		: ProjectiveCamera(pFilm, position, cz, cx, cy, fovy, dis) 
	{

	}

	PinHoleCamera(std::shared_ptr<Film> pFilm, const Vec3& position, const Vec3& lookAt, const Vec3& up, real fovy, real dis)
		: ProjectiveCamera(pFilm, position, lookAt, up, fovy, dis)
	{

	}

	Ray GenerateRay(int pixelX, int pixelY, const Vec2 &u, real offset) const override 
	{

		Vec3 sampledPos = RasterToWorld(Vec3(pixelX + u.x, pixelY + u.y, -1));

		Vec3 dir = (sampledPos - mPos).Norm();
		return Ray(mPos + dir * offset, dir, Inf, 0);
	}

	Vec3 We(const Ray& ray) const override
	{
		real pdfA = 1.0; // for pinhole camera
		real area = mpFilm->mArea;
		real cosTheta = mCz.Dot(ray.mDir);
		real cos2Theta = cosTheta * cosTheta;
		real value = mFilmDistance * mFilmDistance * pdfA / (area * cos2Theta * cos2Theta);
		return Vec3(value, value, value);
	}

	Vec3 Sample_Wi(const Intersection &isect, real *pdfW, Vec3 *wi, Vec3 u) const override 
	{
		*wi = (mPos - isect.mPos);
		real distance = wi->Length();
		*wi = wi->Norm();
		real cosTheta = mCz.Dot(-1 * (*wi));
		*pdfW = 1.f * (distance * distance) / cosTheta;
		return We(Ray(isect.mPos, -1 * (*wi)));
	}


	GYT_FORCE_INLINE real PdfPos() const override 
	{
		return 1.0;
	}

	real PdfDir(const Ray& cameraRay) const override 
	{
		real filmArea = mpFilm->Area();
		real cosTheta = std::abs(mCz.Dot(cameraRay.mDir));
		real cos2Theta = cosTheta * cosTheta;
		return mFilmDistance * mFilmDistance / (filmArea * cos2Theta * cosTheta);
	}

};

GYT_NAMESPACE_END