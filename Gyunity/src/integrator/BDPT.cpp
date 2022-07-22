#include "BDPT.h"
#include "PathVertex.h"
#include "BidirectionalRenderer.h"
#include <mutex>

GYT_NAMESPACE_BEGIN


void BDPT::Render(const Scene& scene, const Camera& camera)
{
	int64 totalPhotons = mSpp * (int64)(camera.GetFilm()->mResX) * (int64)(camera.GetFilm()->mResY);
	int resX = camera.GetFilm()->mResX;
	int resY = camera.GetFilm()->mResY;
	std::atomic<int64> workDone = 0;
	for (int y = 0; y < resY; ++y) {
		std::shared_ptr<Sampler> clonedSampler = mpSampler->Clone(y);
		MemoryPool arena;

		std::vector<PathVertex> lightPath(mMaxDepth + 1);
		std::vector<PathVertex> cameraPath(mMaxDepth + 2);

		for (int x = 0; x < resX; ++x) {
			for (int s = 0; s < mSpp; s++) {
				RandomStateSequence rand(clonedSampler, resX * resY * s + y * resX + x);
				real u = rand();
				real v = rand();
				Ray cameraRay = camera.GenerateRay(x, y, Vec2(u, v), Eps);
				int nLightVertices = BidirectionalRenderer::GenerateLightPath(scene, rand, arena, lightPath, mMaxDepth + 1);

				Vec2 raster;
				Vec3 L = ConnectToCamera(scene, camera, rand, lightPath, nLightVertices, &raster);
				camera.GetFilm()->AddSplat(raster.x, raster.y, L);
				workDone++;
			}
		}
		arena.Reset();
		real percentage = 100.f * workDone / totalPhotons;
		fprintf(stderr, "\rPercentage: %5.2f%%", percentage);
	}
}




Vec3 BDPT::ConnectToCamera(
	const Scene						&scene, 
	const Camera					&camera, 
	StateSequence					&rand,
	const std::vector<PathVertex>	&path, 
	int								s, 
	Vec2							*pRaster)
{
	Vec3 L(0, 0, 0);
	real resX = camera.GetFilm()->mResX;
	real resY = camera.GetFilm()->mResY;
	PathVertex sampled;
	Transform WorldTorRaster = (camera.CameraToRaster * camera.WorldToCamera);

	for (int i = 1; i <= s; ++i) {
		Vec3 raster = WorldTorRaster(path[i].mIsect.mPos);
		bool inScreen = false;
		if (0 <= raster.x && raster.x <= resX && 0 <= raster.y && raster.y <= resY) {
			inScreen = true;
			pRaster->x = raster.x;
			pRaster->y = raster.y;
		}
		if (!BidirectionalRenderer::IsConnectable(scene, path[0].mIsect.mPos, camera.mPos)) {
			inScreen = false;
			L = Vec3(0, 0, 0);
		}
		if (inScreen) {
			real pdfW;
			Vec3 wi;
			Vec3 We = camera.SampleWi(path[i].mIsect, &pdfW, &wi, Vec3(rand(), rand(), rand()));
			Vec3 f = path[i].mIsect.mpBSDF->Evaluate(path[i].mIsect.mOutDir, wi);
			sampled.mIsect.mPos = camera.mPos;
			sampled.mIsect.mNormal = camera.mCz;
			sampled.mThroughput = We / pdfW;
			sampled.mPdfFwd = camera.PdfPos();
			L = We * path[i].mThroughput * f * std::abs(wi.Dot(path[i].mIsect.mNormal)) / pdfW;
		}
	}
	return L;
}

Vec3 BDPT::ConnectLightToCamera(
	const Scene			&scene, 
	const Camera		&camera, 
	const PathVertex	&lightVertex, 
	const PathVertex	&cameraVertex, 
	Vec2				*pRaster)
{
	Vec3 L(0, 0, 0);
	Transform WorldTorRaster = (camera.CameraToRaster * camera.WorldToCamera);
	Vec3 raster = WorldTorRaster(lightVertex.mIsect.mPos);
	bool inScareen = false;
	int resX = camera.GetFilm()->mResX;
	int resY = camera.GetFilm()->mResY;
	if (0 <= raster.x && raster.x <= resX && 0 <= raster.y && raster.y <= resY) {
		inScareen = true;
		pRaster->x = raster.x;
		pRaster->y = raster.y;
	}
	if (BidirectionalRenderer::IsConnectable(scene, lightVertex.mIsect.mPos, cameraVertex.mIsect.mPos)) {
		L = cameraVertex.mThroughput * lightVertex.mThroughput;
	}
	return L;
}

GYT_NAMESPACE_END
