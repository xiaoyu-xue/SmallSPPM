
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
	int photonsPerTask = 10000;
	int remainingPhotons = camera.GetFilm()->mResX * camera.GetFilm()->mResY;
	while(remainingPhotons > 0) {
		if (remainingPhotons > photonsPerTask) {
			tasks.push_back(Task(photonsPerTask));
			remainingPhotons -= photonsPerTask;
		}
		else {
			tasks.push_back(Task(remainingPhotons));
			remainingPhotons -= remainingPhotons;
		}
	}
	std::atomic<int64> workDone = 0;
	int totalWorks = tasks.size() * spp;



	//int phontons = camera.GetFilm()->mResX * camera.GetFilm()->mResY;
	//std::vector<MemoryPool> memoryArenas(ThreadsNumber());
	//std::cout << "Core number: " << ThreadsNumber() << std::endl;
	//for (int i = 0; i < spp; ++i) {
	//	ParallelFor(0, phontons, [&](int p) {
	//		MemoryPool& arena = memoryArenas[ThreadIndex()];
	//		int nVertices = GenerateLightPath(scene, camera, rand, maxDepth, arena);
	//		for (int s = 0; s < nVertices; ++s) {
	//			Vec3 pRaster;
	//			bool inScreen;
	//			Vec3 L = ConnectToCamera(mLightPath[s], s, scene, camera, rand, &pRaster, &inScreen);
	//			if (inScreen) {
	//				L = L * camera.GetFilm()->mResX * camera.GetFilm()->mResY / totalPhotons;
	//				//std::cout << L << std::endl;
	//				camera.GetFilm()->AddSplat(pRaster.x, pRaster.y, L);
	//				//camera.GetFilm()->AddSample(pRaster.x, pRaster.y, L);
	//				arena.Reset();
	//			}
	//		}
	//		workDone += 1;
	//		real percentage = 100.f * workDone / totalPhotons;
	//		fprintf(stderr, "\rPercentage: %5.2f%%", percentage);
	//	});
	//}


	//std::vector<MemoryPool> memoryArenas(ThreadsNumber());
	//std::cout << "Core number: " << ThreadsNumber() << std::endl;
	//ParallelFor(0, totalWorks, [&](int p) {
	//	MemoryPool& arena = memoryArenas[ThreadIndex()];
	//	for (int i = 0; i < tasks[p].mPhotons; ++i) {
	//		int nVertices = GenerateLightPath(scene, camera, rand, maxDepth, arena);
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
	//	}
	//	arena.Reset();
	//	workDone += 1;
	//	real percentage = 100.f * workDone / totalWorks;
	//	fprintf(stderr, "\rPercentage: %5.2f%%", percentage);
	//});

	MemoryPool arena;
	for(int p = 0; p < totalPhotons; ++p) {
		int nVertices = GenerateLightPath(scene, camera, rand, maxDepth, arena);
		for (int s = 0; s < nVertices; ++s) {

			Vec3 pRaster;
			bool inScreen;
			Vec3 L = ConnectToCamera(mLightPath[s], s, scene, camera, rand, &pRaster, &inScreen);
			if (inScreen) {
				L = L * camera.GetFilm()->mResX * camera.GetFilm()->mResY / totalPhotons;
				camera.GetFilm()->AddSplat(pRaster.x, pRaster.y, L);
				arena.Reset();
			}

		}
		real percentage = 100.f * p / totalPhotons;
		fprintf(stderr, "\rPercentage: %5.2f%%", percentage);
	}
}

int LightTracing::GenerateLightPath(const Scene &scene, const Camera &camera, StateSequence& rand, int maxDepth, MemoryPool& arena)
{
	Vec3 dir;
	real lightPdf;
	real pdfDir;
	real pdfA;
	Intersection isect;
	Light *light = scene.SampleOneLight(&lightPdf, rand());
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
	mLightPath[0].mIsect.mPos = isect.mPos;
	mLightPath[0].mIsect.mNormal = isect.mNormal;
	mLightPath[0].mThroughput = Le; // / Pdfdir / Pdfpos * CosTheta;
	mLightPath[0].mPdfFwd = pdfA;
	mLightPath[0].mIsect.mIsDelta = light->IsDeltaLight();
	Vec3 throughput = mLightPath[0].mThroughput * cosTheta / lightPdf / pdfA / pdfDir;
	Ray ray(isect.mPos + dir * RayEps, dir);

	int nVetices = Trace(ray, throughput, pdfDir, rand, 1, maxDepth, scene, camera, arena);

	return nVetices;
}

int LightTracing::Trace(const Ray& ray, Vec3 throughput, real pdfFwd, StateSequence& rand, int depth, int maxDepth, const Scene& scene, const Camera& camera, MemoryPool& arena)
{
	real pdfdir = pdfFwd;
	Ray r = ray;
	int bound = depth;
	real pdfW = pdfFwd;

	while (1) {
		Intersection& isect = mLightPath[bound].mIsect;
		if (!scene.Intersect(r, &isect)) break;

		scene.QueryIntersectionInfo(r, &isect);
		isect.ComputeScatteringFunction(arena, TransportMode::Importance);

		mLightPath[bound].mThroughput = throughput;
		mLightPath[bound].mPdfFwd = pdfW;

		++bound;
		if (bound >= maxDepth + 1) break;

		Vec3 wo;
		Vec3 f = isect.mpBSDF->Sample(-1 * r.mDir, &wo, &pdfW, Vec3(rand(), rand(), rand()));
		Normalize(wo);
		throughput = throughput * f * (wo.Dot(isect.mNormal)) / pdfW;

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
	
	if (vertex.mIsect.mIsDelta) return Vec3(0.0, 0.0, 0.0);

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

	Ray ray(vertex.mIsect.mPos, hitPointToCam, 0.f, ray_tmax);
	Intersection isc;

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
	Vec3 L = We * vertex.mThroughput * f * wi.Dot(vertex.mIsect.mNormal) / pdfW;

	return L;
}


GYT_NAMESPACE_END
