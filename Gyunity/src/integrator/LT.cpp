
#include "LT.h"
#include "PathVertex.h"

GYT_NAMESPACE_BEGIN


struct Task
{
	int mPhotons;
	Task(int n) : mPhotons(n) {}
};

void LightTracing::Render(const Scene& scene, const Camera& camera)
{
	RandomStateSequence rand(mpSampler, 123);
	int totalPhotons = spp * camera.GetFilm()->mResX * camera.GetFilm()->mResY;
	std::cout << totalPhotons << std::endl;

	std::vector<Task> tasks;
	int photonsPerTask = 100000;
	int remainingPhotons = camera.GetFilm()->mResX * camera.GetFilm()->mResY;
	for (int s = 0; s < spp; ++s) {
		while (remainingPhotons > 0) {
			if (remainingPhotons > photonsPerTask) {
				tasks.push_back(Task(photonsPerTask));
				remainingPhotons -= photonsPerTask;
			}
			else {
				tasks.push_back(Task(remainingPhotons));
				remainingPhotons -= remainingPhotons;
			}
		}
	}
	std::atomic<int64> workDone = 0;
	int totalWorks = tasks.size();



	//int phontons = camera.GetFilm()->mResX * camera.GetFilm()->mResY;
	//std::vector<MemoryPool> memoryArenas(ThreadsNumber());
	//std::cout << "Core number: " << ThreadsNumber() << std::endl;
	//for (int i = 0; i < spp; ++i) {
	//	ParallelFor(0, phontons, [&](int p) {
	//		MemoryPool& arena = memoryArenas[ThreadIndex()];
	//		int nVertices = GenerateLightPath(scene, camera, rand, arena);
	//		for (int s = 0; s < nVertices; ++s) {
	//			Vec3 pRaster;
	//			bool inScreen;
	//			Vec3 L = ConnectToCamera(mLightPath[s], s, scene, camera, rand, &pRaster, &inScreen);
	//			if (inScreen) {
	//				L = L * camera.GetFilm()->mResX * camera.GetFilm()->mResY / totalPhotons;
	//				//std::cout << L << std::endl;
	//				camera.GetFilm()->AddSplat(pRaster.x, pRaster.y, L);
	//				//camera.GetFilm()->AddSample(pRaster.x, pRaster.y, L);
	//				
	//			}
	//		}
	//		arena.Reset();
	//		workDone += 1;
	//		real percentage = 100.f * workDone / totalPhotons;
	//		fprintf(stderr, "\rPercentage: %5.2f%%", percentage);
	//	});
	//}


	//std::vector<MemoryPool> memoryArenas(ThreadsNumber());
	//GYT_Print("Core number: {}\n", ThreadsNumber());
	//GYT_Print("Total works: {}\n", totalWorks);
	//ParallelFor(0, totalWorks, [&](int p) {
	//	MemoryPool& arena = memoryArenas[ThreadIndex()];
	//	for (int i = 0; i < tasks[p].mPhotons; ++i) {
	//		std::vector<PathVertex> lightPath(maxDepth + 1);
	//		int nVertices = GenerateLightPath(scene, camera, rand, arena, lightPath);
	//		for (int s = 0; s < nVertices; ++s) {
	//			Vec3 pRaster;
	//			bool inScreen;
	//			Vec3 L = ConnectToCamera(lightPath[s], s, scene, camera, rand, &pRaster, &inScreen);
	//			if (inScreen) {
	//				//L = L * camera.GetFilm()->mResX * camera.GetFilm()->mResY / totalPhotons;
	//				L = L / spp;
	//				//std::cout << L << std::endl;
	//				camera.GetFilm()->AddSplat(pRaster.x, pRaster.y, L);
	//				//camera.GetFilm()->AddSample(pRaster.x, pRaster.y, L);
	//				
	//			}
	//		}
	//	}
	//	arena.Reset();
	//	workDone += 1;
	//	real percentage = 100.f * workDone / totalWorks;
	//	fprintf(stderr, "\rPercentage: %5.2f%%", percentage);
	//});

	MemoryPool arena;
	for (int p = 0; p < totalPhotons; ++p) {
		std::vector<PathVertex> lightPath(maxDepth + 1);
		int nVertices = GenerateLightPath(scene, camera, rand, arena, lightPath);
		for (int s = 0; s < nVertices; ++s) {
			Vec3 pRaster;
			bool inScreen;
			Vec3 L = ConnectToCamera(lightPath[s], s, scene, camera, rand, &pRaster, &inScreen);
			if (inScreen) {
				L = L * camera.GetFilm()->mResX * camera.GetFilm()->mResY / totalPhotons;
				camera.GetFilm()->AddSplat(pRaster.x, pRaster.y, L);
			}
		}
		arena.Reset();
		real percentage = 100.f * p / totalPhotons;
		fprintf(stderr, "\rPercentage: %5.2f%%", percentage);
	}
}

int LightTracing::GenerateLightPath(const Scene &scene, const Camera &camera, StateSequence& rand, MemoryPool& arena, std::vector<PathVertex> &lightPath)
{
	Vec3 dir;
	real pdfLight;
	real pdfDir;
	real pdfA;
	Intersection isect;
	Light *light = scene.SampleOneLight(&pdfLight, rand());
	Vec3 Le;
	if (light->IsAreaLight()) {
		Le = light->Emission();
		light->SampleOnLight(&isect, &dir, &pdfA, &pdfDir, Vec2(rand(), rand()), Vec2(rand(), rand()));
	}
	else {
		Le = light->Emission();
		light->SampleLight(&isect, &dir, &pdfA, &pdfDir, Vec2(rand(), rand()), Vec2(rand(), rand()));
	}
	real cosTheta = isect.mNormal.Dot(dir);
	lightPath[0].mIsect.mPos = isect.mPos;
	lightPath[0].mIsect.mNormal = isect.mNormal;
	lightPath[0].mThroughput = Le; // / Pdfdir / Pdfpos * CosTheta;
	lightPath[0].mPdfFwd = pdfA;
	lightPath[0].mIsDelta = light->IsDeltaLight();

	//GYT_Print("cosTheta: {}, pdfLight: {}, pdfA: {}, pdfDir: {}\n", cosTheta, pdfLight, pdfA, pdfDir);

	Vec3 throughput = lightPath[0].mThroughput * cosTheta / pdfLight / pdfA / pdfDir;
	Ray ray(isect.mPos + dir * RayEps, dir);

	int nVetices = Trace(ray, throughput, pdfDir, rand, 1, scene, camera, arena, lightPath);

	return nVetices;
}

int LightTracing::Trace(const Ray& ray, Vec3 throughput, real pdfFwd, StateSequence& rand, int depth, const Scene& scene, const Camera& camera, MemoryPool& arena, std::vector<PathVertex>& lightPath)
{
	real pdfdir = pdfFwd;
	Ray r = ray;
	int bound = depth;
	real pdfW = pdfFwd;

	while (1) {
		Intersection& isect = lightPath[bound].mIsect;
		if (!scene.Intersect(r, &isect)) break;

		scene.QueryIntersectionInfo(r, &isect);
		isect.ComputeScatteringFunction(arena, TransportMode::Importance);

		lightPath[bound].mIsDelta = isect.mpBSDF->IsDelta();
		lightPath[bound].mThroughput = throughput;
		lightPath[bound].mPdfFwd = pdfW;

		GYT_Print("{}\n", throughput);

		++bound;
		if (bound >= maxDepth + 1) break;

		Vec3 wo;
		Vec3 f = isect.mpBSDF->Sample(-1 * r.mDir, &wo, &pdfW, Vec3(rand(), rand(), rand()));
		Normalize(wo);
		throughput = throughput * f * std::abs(wo.Dot(isect.mNormal)) / pdfW;

		r.mOrig = isect.mPos + wo * RayEps;
		r.mDir = wo;
	}

	return bound - 1;
}



Vec3 LightTracing::WorldToScreen(const Camera & camera, const Vec3& vertex, bool* isInScreen) const
{
	*isInScreen = false;
	Transform WorldToRaster = camera.CameraToRaster * camera.WorldToCamera;
	Vec3 rasterPos = WorldToRaster(vertex);
	if (rasterPos.x >= 0 && rasterPos.x <= camera.GetFilm()->mResX &&
		rasterPos.y >= 0 && rasterPos.y <= camera.GetFilm()->mResY) {
		*isInScreen = true;
	}
	return rasterPos;
}


Vec3 LightTracing::ConnectToCamera(const PathVertex& vertex, int s, const Scene& scene, const Camera& camera, StateSequence& rand, Vec3* pRaster, bool* inScreen) {
	
	if (vertex.mIsDelta) 
		return Vec3(0.0, 0.0, 0.0);

	*inScreen = true;
	Vec3 hitPointToCam = camera.mPos - vertex.mIsect.mPos;
	real ray_tmax = hitPointToCam.Length();
	Normalize(hitPointToCam);
	if (camera.mCz.Dot(-1 * hitPointToCam) < 0) {
		*inScreen = false;
		return Vec3(0.0, 0.0, 0.0);
	}

	bool isInScreen;
	*pRaster = WorldToScreen(camera, vertex.mIsect.mPos, &isInScreen);

	if (!isInScreen) {
		*inScreen = false;
		return Vec3(0.0, 0.0, 0.0);
	}

	Ray ray(vertex.mIsect.mPos + hitPointToCam * RayEps, hitPointToCam, 0.f, ray_tmax);

	if (scene.Intersect(ray)) {
		* inScreen = false;
		return Vec3(0.0, 0.0, 0.0);
	}
	if (s == 0) {
		return vertex.mThroughput; //see the light source directly
	}
	real pdfW;
	Vec3 wi;
	Vec3 We = camera.Sample_Wi(vertex.mIsect, &pdfW, &wi, Vec3(rand(), rand(), rand()));
	Vec3 f = vertex.mIsect.mpBSDF->Evaluate(vertex.mIsect.mOutDir, wi);
	Vec3 L = We * vertex.mThroughput * f * std::abs(wi.Dot(vertex.mIsect.mNormal)) / pdfW;

	return L;
}


GYT_NAMESPACE_END
