#include "Camera.h"
#include "visual/Film.h"

GYT_NAMESPACE_BEGIN

ProjectiveCamera::ProjectiveCamera(const std::shared_ptr<Film>& pFilm, const Vec3& position,
	const Vec3& pCz, const Vec3& pCx, const Vec3& pCy, real pFovy, real dis) :
	Camera(pFilm), pos(position), cx(pCx), cy(pCy), cz(pCz), fovy(pFovy), filmDistance(dis) {

	Initialize();

	CameraToWorld =
		Transform(cx.x, cy.x, cz.x, pos.x,
			cx.y, cy.y, cz.y, pos.y,
			cx.z, cy.z, cz.z, pos.z,
			0, 0, 0, 1.f);

	WorldToCamera = Inverse(CameraToWorld);

	CameraToNDC = Transform::Perspective(fovy, film->mAspect, filmDistance, filmDistance, 100000.f);

	NDCToCamera = Inverse(CameraToNDC);

	NDCToRaster =
		Transform::Scale(real(film->mResX), real(film->mResY), 1) *
		Transform::Scale(0.5f, -0.5f, 1) *
		Transform::Translate(Vec3(1, -1, 0));

	RasterToNDC = Inverse(NDCToRaster);

	CameraToRaster = NDCToRaster * CameraToNDC;

	RasterToCamera = NDCToCamera * RasterToNDC;

	RasterToWorld = CameraToWorld * RasterToCamera;

}

void ProjectiveCamera::Initialize() {
	Vec3 filmCenter = pos + czz * filmDistance;
	//std::cout << "film cetner: " << filmCenter << std::endl;
	film->mHeight = filmDistance * std::tan(fovy * 0.5f * PI / 180) * 2.f;
	film->mWidth = film->mHeight * film->mAspect;
	film->mArea = film->mWidth * film->mHeight;

	film->mLU = filmCenter + cy * film->mHeight * 0.5 - cx * film->mWidth * 0.5;
	film->mLL = filmCenter - cy * film->mHeight * 0.5 - cx * film->mWidth * 0.5;
	film->mRU = filmCenter + cy * film->mHeight * 0.5 + cx * film->mWidth * 0.5;
	film->mRL = filmCenter - cy * film->mHeight * 0.5 + cx * film->mWidth * 0.5;

	GYT_Print("LL: {} \n", film->mLL);
	GYT_Print("LL: {} \n", film->mLU);
	GYT_Print("LL: {} \n", film->mRL);
	GYT_Print("LL: {} \n", film->mRU);

}


GYT_NAMESPACE_END