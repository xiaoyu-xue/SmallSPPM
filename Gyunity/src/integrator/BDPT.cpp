#include "BDPT.h"
#include "PathVertex.h"
#include "BidirectionalRenderer.h"
#include "system/Threading.h"
#include "sampler/SamplerEnum.h"
#include <mutex>

GYT_NAMESPACE_BEGIN

void BDPT::Render(const Scene& scene, const Camera& camera)
{
	int64 resX = camera.GetFilm()->mResX;
	int64 resY = camera.GetFilm()->mResY;
	int64 totalPhotons = mSpp * (int64)(camera.GetFilm()->mResX) * (int64)(camera.GetFilm()->mResY);
	std::vector<MemoryPool> memoryArenas(ThreadsNumber());
	std::atomic<int64> workDone = 0;

	ParallelFor(int64(0), resY, [&](int64 y) {

		std::shared_ptr<Sampler> clonedSampler = mpSampler->Clone(y);
		MemoryPool& arena = memoryArenas[ThreadIndex()];

		std::vector<PathVertex> lightPath(mMaxDepth + 1);
		std::vector<PathVertex> cameraPath(mMaxDepth + 2);
		for (int64 x = 0; x < resX; ++x) {
			for (int s = 0; s < mSpp; ++s) {
				RandomStateSequence rand(clonedSampler, resX * resY * s + y * resX + x);

				Ray cameraRay = camera.GenerateRay(x, y, Vec2(rand(), rand()));

				//int nLightVertices = BidirectionalRenderer::GenerateLightPath(scene, rand, arena, lightPath, mMaxDepth);
				int nCameraVertices = BidirectionalRenderer::GenerateCameraPath(scene, camera, rand, arena, cameraPath, cameraRay, mMaxDepth);

				real u = mpSamplerEnum->SampleX(x, rand());
				real v = mpSamplerEnum->SampleY(y, rand());
				//for (int s = 0; s < nLightVertices; ++s) {
				//	Vec3 pRaster;
				//	bool inScreen;
				//	Vec3 L = ConnectToCamera(scene, camera, rand, lightPath[s], s, &pRaster, &inScreen);
				//	if (inScreen) {
				//		//L = L / mSpp;
				//		L = L * camera.GetFilm()->mResX * camera.GetFilm()->mResY / totalPhotons;
				//		camera.GetFilm()->AddSplat(pRaster.x, pRaster.y, L);
				//	}
				//}
				for (int t = 1; t <= 1; ++t) {
					Vec3 L = ConnectToLight(scene, rand, cameraPath, t);
					camera.GetFilm()->AddSample(x + u, y + v, L);
				}
				workDone++;
			}
		}
		arena.Reset();
		real percentage = 100.f * workDone / totalPhotons;
		fprintf(stderr, "\rPercentage: %5.2f%%", percentage);
	});
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

Vec3 BDPT::ConnectToCamera(const Scene& scene, const Camera& camera, StateSequence& rand, const PathVertex& vertex, int s, Vec3* pRaster, bool* inScreen) const
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
	Intersection isect;
	real pdf;
	Vec3 L(0, 0, 0);
	pLight = scene.SampleOneLight(&pdfLight, rand());
	pLight->SampleOnePoint(&isect, &pdf, Vec2(rand(), rand()));

	Vec3 wi = (isect.mPos - path[t].mIsect.mPos);
	Ray shadowRay(path[t].mIsect.mPos + RayEps, wi.Length(), 0);
	wi.Normalize();
	if (!scene.Intersect(shadowRay)) {
		if (t == 1) {
			return pLight->Emission();
		}
	}
	Vec3 f = path[t].mIsect.mpBSDF->Evaluate(path[t].mIsect.mOutDir, wi);
	//L = path[t].mThroughput * f * std::abs(wi.Dot(path[t].mIsect.mNormal)) * pLight->Emission() / pdfLight / pdf;
	//Vec3 L = path[t].mThroughput * DirectIllumination(scene, path[t].mIsect, rand(), Vec2(rand(), rand()), Vec3(rand(), rand(), rand()), rand) / pdfLight / pdf ;
	return L;
}


GYT_NAMESPACE_END
