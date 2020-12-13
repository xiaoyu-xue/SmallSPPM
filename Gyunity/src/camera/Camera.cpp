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

void ProjectiveCamera::Initialize() {
	Vec3 filmCenter = pos + czz * filmDistance;
	//std::cout << "film cetner: " << filmCenter << std::endl;
	film->height = filmDistance * std::tan(fovy * 0.5f * PI / 180) * 2.f;
	film->width = film->height * film->aspect;
	film->area = film->width * film->height;

	film->LU = filmCenter + cy * film->height * 0.5 - cx * film->width * 0.5;
	film->LL = filmCenter - cy * film->height * 0.5 - cx * film->width * 0.5;
	film->RU = filmCenter + cy * film->height * 0.5 + cx * film->width * 0.5;
	film->RL = filmCenter - cy * film->height * 0.5 + cx * film->width * 0.5;

	GYT_Print("LL: {} \n", film->LL);
	GYT_Print("LL: {} \n", film->LU);
	GYT_Print("LL: {} \n", film->RL);
	GYT_Print("LL: {} \n", film->RU);

}


GYT_NAMESPACE_END