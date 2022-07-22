#pragma once
#include <vector>
#include "ForwardDecl.h"
#include "Integrator.h"
#include "sampler/Sampler.h"
#include "PathVertex.h"
#include <mutex>

GYT_NAMESPACE_BEGIN

class LightTracing : public Integrator {

public:
	LightTracing(std::shared_ptr<Sampler> pSampler, int maxDepth, int spp) :
		mpSampler(pSampler), mSpp(spp), mMaxDepth(maxDepth) { }

	void Render(const Scene& scene, const Camera& camera) override {

		int64 totalPhotons = mSpp * (int64)(camera.GetFilm()->mResX) * (int64)(camera.GetFilm()->mResY);
		int64 photonsPerRenderStage = (int64)(camera.GetFilm()->mResX) * (int64)(camera.GetFilm()->mResY);
		std::atomic<int64> workDone = 0;

		int64 resX = camera.GetFilm()->mResX;
		int64 resY = camera.GetFilm()->mResY;
		std::vector<MemoryPool> memoryArenas(ThreadsNumber());
		ParallelFor(int64(0), resY, [&](int64 y) {
			std::shared_ptr<Sampler> clonedSampler = mpSampler->Clone(y);
			MemoryPool& arena = memoryArenas[ThreadIndex()];
			for (int64 x = 0; x < resX; ++x) {
				for (int p = 0; p < mSpp; ++p) {
					RandomStateSequence rand(clonedSampler, resX * resY * p + y * resX + x);
					std::vector<PathVertex> lightPath(mMaxDepth + 1);
					int nVertices = GenerateLightPath(scene, camera, rand, arena, lightPath);
					for (int s = 0; s < nVertices; ++s) {
						Vec3 pRaster;
						bool inScreen;
						Vec3 L = ConnectToCamera(scene, camera, rand, lightPath[s], s, &pRaster, &inScreen);
						if (inScreen) {
							//L = L / mSpp;
							L = L * camera.GetFilm()->mResX * camera.GetFilm()->mResY / totalPhotons;
							camera.GetFilm()->AddSplat(pRaster.x, pRaster.y, L);
						}
					}
					workDone++;
				}
			}
			arena.Reset();
			real percentage = 100.f * workDone / totalPhotons;
			fprintf(stderr, "\rPercentage: %5.2f%%", percentage);
		});
	}

	int Trace(const Scene& scene, const Camera& camera, const Ray& r, MemoryPool& arena, StateSequence& rand,
		int depth, Vec3 throughput, real pdfFwd, std::vector<PathVertex>& lightPath) const {
		
		Ray ray = r;
		int bounce = depth;
		real pdfDir = pdfFwd;
		while (true) {
			Intersection& isect = lightPath[bounce].mIsect;

			if (!scene.Intersect(ray, &isect)) break;

			scene.QueryIntersectionInfo(ray, &isect);
			isect.ComputeScatteringFunction(arena, TransportMode::Importance);

			lightPath[bounce].mIsDelta = isect.mpBSDF->IsDelta();
			lightPath[bounce].mPdfFwd = pdfDir;
			lightPath[bounce].mThroughput = throughput;

			bounce++;
			if (bounce >= mMaxDepth) break;

			Vec3 wi;
			Vec3 f = isect.mpBSDF->Sample(-1 * ray.mDir, &wi, &pdfDir, Vec3(rand(), rand(), rand()));
			ray = Ray(isect.mPos + wi * RayEps, wi);

			throughput = throughput * f * std::abs(Dot(wi, isect.mNormal)) / pdfDir;

		}
		return bounce;
	}

	int GenerateLightPath(const Scene& scene, const Camera& camera, StateSequence& rand, MemoryPool& arena, 
		std::vector<PathVertex>& lightPath) const {
		
		real pdfLight;
		real pdfA;
		real pdfDir;
		Vec3 Le;
		Vec3 dir;
		Intersection isect;
		Light* light = scene.SampleOneLight(&pdfLight, rand());
		if (light->IsAreaLight()) {
			Le = light->Emission();
			light->SampleOnLight(&isect, &dir, &pdfA, &pdfDir, Vec2(rand(), rand()), Vec2(rand(), rand()));
		}
		else {
			Le = light->Emission();
			light->SampleLight(&isect, &dir, &pdfA, &pdfDir, Vec2(rand(), rand()), Vec2(rand(), rand()));
		}
		lightPath[0].mIsect = isect;
		lightPath[0].mPdfFwd = pdfA;
		lightPath[0].mThroughput = Le;
		lightPath[0].mIsDelta = light->IsDeltaLight();
		real cosTheta = std::abs(Dot(isect.mNormal, dir));
		Vec3 throughput = Le * cosTheta / pdfLight / pdfA / pdfDir;
		Ray ray(isect.mPos + dir * RayEps, dir);

		int nVertices = Trace(scene, camera, ray, arena, rand, 1, throughput, pdfDir, lightPath);
		


		return nVertices;
	}

	Vec3 WorldToScreen(const Camera& camera, const Vec3& point, bool* inScreen) const {
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


	Vec3 ConnectToCamera(const Scene& scene, const Camera& camera, StateSequence& rand,
		const PathVertex& vertex, int s, Vec3* pRaster, bool* inScreen) const {
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

private:
	int mSpp, mMaxDepth;
	std::shared_ptr<Sampler> mpSampler;

	//debug
	mutable std::mutex mMutex;
};

GYT_NAMESPACE_END