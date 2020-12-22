#include "Camera.h"
#include "visual/Film.h"

GYT_NAMESPACE_BEGIN

void ProjectiveCamera::Initialize() {

	mCz = -(mPos - mLookAt).Norm(); // -1 * dir in raytracer
	mCx = -Cross(mUp, mCz).Norm();
	mCy = -Cross(mCz, mCx);

	CameraToWorld =
		Transform(mCx.x, mCy.x, mCz.x, mPos.x,
			mCx.y, mCy.y, mCz.y, mPos.y,
			mCx.z, mCy.z, mCz.z, mPos.z,
			0, 0, 0, 1.f);

	//Should get the same result
	//CameraToWorld = Transform::LookAt(mPos, mLookAt, mUp);

	WorldToCamera = Inverse(CameraToWorld);

	CameraToNDC = Transform::Perspective(mFovy, mpFilm->mAspectRatio, mFilmDistance, mFilmDistance, 100000.f);

	NDCToCamera = Inverse(CameraToNDC);

	NDCToRaster =
		Transform::Scale(real(mpFilm->mResX), real(mpFilm->mResY), 1) *
		Transform::Scale(0.5f, -0.5f, 1) *
		Transform::Translate(Vec3(1, -1, 0));

	RasterToNDC = Inverse(NDCToRaster);

	CameraToRaster = NDCToRaster * CameraToNDC;

	RasterToCamera = NDCToCamera * RasterToNDC;

	RasterToWorld = CameraToWorld * RasterToCamera;


	Vec3 filmCenter = mPos + mCz * mFilmDistance;
	GYT_Print("Film center: {}\n", filmCenter);
	mpFilm->mHeight = mFilmDistance * std::tan(mFovy * 0.5f * PI / 180) * 2.f;
	mpFilm->mWidth = mpFilm->mHeight * mpFilm->mAspectRatio;
	mpFilm->mArea = mpFilm->mWidth * mpFilm->mHeight;

	mpFilm->mLU = filmCenter + mCy * mpFilm->mHeight * 0.5 - mCx * mpFilm->mWidth * 0.5;
	mpFilm->mLL = filmCenter - mCy * mpFilm->mHeight * 0.5 - mCx * mpFilm->mWidth * 0.5;
	mpFilm->mRU = filmCenter + mCy * mpFilm->mHeight * 0.5 + mCx * mpFilm->mWidth * 0.5;
	mpFilm->mRL = filmCenter - mCy * mpFilm->mHeight * 0.5 + mCx * mpFilm->mWidth * 0.5;

	GYT_Print("LL: {} \n", mpFilm->mLL);
	GYT_Print("LL: {} \n", mpFilm->mLU);
	GYT_Print("LL: {} \n", mpFilm->mRL);
	GYT_Print("LL: {} \n", mpFilm->mRU);

}


ProjectiveCamera::ProjectiveCamera(std::shared_ptr<Film> pFilm, const Vec3& position, const Vec3& cz, const Vec3& cx, const Vec3& cy, real pFovy, real dis)
	: Camera(pFilm), mPos(position), mCx(cx), mCy(cy), mCz(cz), mFovy(pFovy), mFilmDistance(dis)
{
	Initialize();
}

ProjectiveCamera::ProjectiveCamera(std::shared_ptr<Film> pFilm, const Vec3& position, const Vec3& lookAt, const Vec3& up, real fovy, real dis)
	: Camera(pFilm), mPos(position), mLookAt(lookAt), mUp(up), mFovy(fovy), mFilmDistance(dis)
{
	Initialize();
}

GYT_NAMESPACE_END