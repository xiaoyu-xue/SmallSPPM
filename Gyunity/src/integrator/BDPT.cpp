#include "BDPT.h"
#include "PathVertex.h"
#include "BidirectionalRenderer.h"
#include "system/Threading.h"
#include <mutex>

GYT_NAMESPACE_BEGIN


void BDPT::Render(const Scene& scene, const Camera& camera)
{
	int64 totalPhotons = mSpp * (int64)(camera.GetFilm()->mResX) * (int64)(camera.GetFilm()->mResY);
	int resX = camera.GetFilm()->mResX;
	int resY = camera.GetFilm()->mResY;
	std::vector<MemoryPool> memoryArenas(ThreadsNumber());
	std::atomic<int64> workDone = 0;
	ParallelFor(int64(0), (int64)resY, [&](int64 y) {
	//for(int y = 0; y < resY; ++y){
		std::shared_ptr<Sampler> clonedSampler = mpSampler->Clone(y);
		MemoryPool& arena = memoryArenas[ThreadIndex()];
		std::vector<PathVertex> lightPath(mMaxDepth + 1);
		std::vector<PathVertex> cameraPath(mMaxDepth + 2);
		for (int x = 0; x < resX; ++x) {
			for (int s = 0; s < mSpp; s++) {
				RandomStateSequence rand(clonedSampler, resX * resY * s + y * resX + x);
				real u = rand();
				real v = rand();
				Ray cameraRay = camera.GenerateRay(x, y, Vec2(u, v), Eps);
				//int nLightVertices = BidirectionalRenderer::GenerateLightPath(scene, rand, arena, lightPath, mMaxDepth + 1);
				int nCameraVertices = BidirectionalRenderer::GenerateCameraPath(camera, rand, arena, cameraPath, cameraRay, mMaxDepth + 2);
				Vec2 raster;
				Vec2 pRaster;
				bool inScreen;

				//for (int s = 0; s < nLightVertices; s++) {
				//	Vec3 L = ConnectToCameraV2(scene, camera, rand, lightPath[s], s, &pRaster, &inScreen);
				//	L = L / mSpp;
				//	camera.GetFilm()->AddSplat(pRaster.x, pRaster.y, L);
				//}
				for (int t = 0; t < nCameraVertices; t++) {
					Vec3 L = ConnectToLight(scene, rand, cameraPath, t);
					L = L / mSpp;
				}
				workDone++;
			}
		}
		arena.Reset();
		real percentage = 100.f * workDone / totalPhotons;
		fprintf(stderr, "\rPercentage: %5.2f%%", percentage);
	});
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
			//sampled.mIsect.mPos = camera.mPos;
			//sampled.mIsect.mNormal = camera.mCz;
			//sampled.mThroughput = We / pdfW;
			//sampled.mPdfFwd = camera.PdfPos();
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

Vec3 BDPT::WorldToScreen(const Camera& camera, const Vec3& point, bool* inScreen) const
{
	Transform World2Raster = camera.CameraToRaster * camera.WorldToCamera;
	Vec3 rasterPoint = World2Raster(point);
	Vec2i rasterPointPrim = Vec2i(int(rasterPoint.x), int(rasterPoint.y));
	if (rasterPointPrim.x >= 0 && rasterPointPrim.x <= camera.GetFilm()->mResX &&
		rasterPointPrim.y >= 0 && rasterPointPrim.y <= camera.GetFilm()->mResY) {
		*inScreen = true;
	}
	else {
		*inScreen = false;
	}
	return Vec3(rasterPointPrim.x, rasterPointPrim.y, 0);
}

Vec3 BDPT::ConnectToCameraV2(const Scene& scene, const Camera& camera, StateSequence& rand, const PathVertex& vertex, int s, Vec3* pRaster, bool* inScreen)
{
	Vec3 pos = vertex.mIsect.mPos;
	Vec3 hitToCameraDir = camera.mPos - pos;
	real hitToCamLength = hitToCameraDir.Length();
	Normalize(hitToCameraDir);
	if (Dot(camera.mCz, -1 * hitToCameraDir) < 0) {
		*inScreen = false;
		return Vec3();
	}

	*pRaster = WorldToScreen(camera, pos, inScreen);

	Ray ray(pos, hitToCameraDir, hitToCamLength);
	if (scene.Intersect(ray))
	{
		*inScreen = false;
	}

	if (s == 0) {
		return vertex.mThroughput;
	}

	real pdfW;
	Vec3 wi;
	Vec3 We = camera.SampleWi(vertex.mIsect, &pdfW, &wi, Vec3(rand(), rand(), rand()));
	Vec3 f = vertex.mIsect.mpBSDF->Evaluate(vertex.mIsect.mOutDir, wi);
	Vec3 L = We * f * vertex.mThroughput * std::abs(Dot(wi, vertex.mIsect.mNormal)) / pdfW;

	return L;
}

Gyunity::Vec3 BDPT::ConnectToLight(const Scene& scene, StateSequence& rand, const std::vector<PathVertex>& path, int t)
{
	Light* pLight;
	real pdfLight;
	pLight = scene.SampleOneLight(&pdfLight, rand());
	real pdfPos, pdfDir;
	Intersection isect;
	Vec3 dir;
	pLight->SampleOnLight(&isect, &dir, &pdfPos, &pdfDir, Vec2(rand(), rand()), Vec2(rand(), rand()));
	Vec3 wi = -1 * dir;
	Vec3 f = path[t].mIsect.mpBSDF->Evaluate(path[t].mIsect.mOutDir, wi);
	Ray shadowRay(path[t].mIsect.mPos, dir.Norm());
	if (scene.Intersect(shadowRay)) {
		return Vec3(0, 0, 0);
	}
	Vec3 L = path[t].mThroughput * f * wi.Dot(path[t].mIsect.mNormal) * pLight->Emission() / pdfDir / pdfLight;
	return L;
}

GYT_NAMESPACE_END
