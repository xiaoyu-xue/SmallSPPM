#include "PathTracingV2.h"
#include "PathVertex.h"
#include "system/Threading.h"
#include "system/Memory.h"
#include "sampler/SamplerEnum.h"

GYT_NAMESPACE_BEGIN

real ConvertSolidToArea(real pdfW, const PathVertex& vertex, const PathVertex& nxt) {
	Vec3 dir = nxt.mIsect.mPos - vertex.mIsect.mPos;
	double dist = dir.Length();
	double dist2 = dist * dist;
	if (dist2 == 0) return 0;
	dir.Normalize();
	double cosTheta = std::abs(dir.Dot(nxt.mIsect.mNormal));
	return pdfW * cosTheta / dist2;
}

void PathTracingV2::Render(const Scene& scene, const Camera& camera)
{
	int64 resX = camera.GetFilm()->mResX;
	int64 resY = camera.GetFilm()->mResY;
	int64 totalPhotons = mSpp * (int64)(camera.GetFilm()->mResX) * (int64)(camera.GetFilm()->mResY);
	std::vector<MemoryPool> memoryArenas(ThreadsNumber());
	std::atomic<int64> workDone = 0;

	ParallelFor(int64(0), resY, [&](int64 y) {
		std::shared_ptr<Sampler> clonedSampler = mpSampler->Clone(y);
		MemoryPool& arena = memoryArenas[ThreadIndex()];
		for (int64 x = 0; x < resX; ++x) {
			std::vector<PathVertex> cameraPath(mMaxDepth + 2);
			for (int spp = 0; spp < mSpp; spp++) {
				RandomStateSequence rand(clonedSampler, resX * resY * spp + y * resX + x);
				real u = mpSamplerEnum->SampleX(x, rand());
				real v = mpSamplerEnum->SampleY(y, rand());
				Vec2 pixelSample(u, v);
				Ray cameraRay = camera.GenerateRay(x, y, pixelSample);
				//int nCameraVertices = ConstructCameraPath(scene, camera, rand, arena, cameraRay, cameraPath, 15);
				int nCameraVertices = GenerateCameraPath(scene, camera, rand, arena, cameraPath, cameraRay, mMaxDepth);
				Vec3 L(0, 0, 0);

				for (int t = 1; t < nCameraVertices; t++) {
					L += ConnectToLight(scene, rand, cameraPath, t);
				}

				//Vec3 L = Li(scene, rand, arena, ray);
				camera.GetFilm()->AddSample(x + u, y + v, L);
				workDone++;
			}
			arena.Reset();
		}
		real percentage = 100.f * workDone / totalPhotons;
		fprintf(stderr, "\rPercentage: %5.2f%%", percentage);
	});
}

int PathTracingV2::GenerateCameraPath(const Scene& scene, const Camera& camera, StateSequence& rand, MemoryPool& arena, std::vector<PathVertex>& cameraPath, const Ray& cameraRay, int maxDepth)
{
	if (maxDepth == 0) return 0;
	Vec3 throughput(1.f, 1.f, 1.f);
	cameraPath[0].mIsect.mPos = camera.mPos;
	cameraPath[0].mIsect.mOutDir = -1 * cameraRay.mDir;
	cameraPath[0].mIsect.mNormal = camera.mCz;
	cameraPath[0].mIsect.mGeometryNormal = camera.mCz;
	cameraPath[0].mPdfFwd = camera.PdfPos();
	cameraPath[0].mThroughput = throughput;
	real pdfDir = camera.PdfDir(cameraRay);
	const Camera* pCam = &camera;
	cameraPath[0].mpCamera = const_cast<Camera*>(pCam);
	int nCameraVertices = Trace(scene, arena, rand, cameraRay, 1, throughput, pdfDir, cameraPath, maxDepth - 1, TransportMode::Radiance);
	return nCameraVertices;
}

int PathTracingV2::Trace(const Scene& scene, MemoryPool& arena, StateSequence& rand, const Ray& r, int depth, Vec3 throughput, real pdfFwd, std::vector<PathVertex>& path, int maxDepth, TransportMode mode)
{
	int bounce = 1;
	Ray ray = r;
	real pdfW = pdfFwd;
	while (true) {
		PathVertex& prev = path[bounce - 1];
		PathVertex& vertex = path[bounce];
		Intersection& isect = path[bounce].mIsect;
		if (!scene.Intersect(ray, &isect)) {
			break;
		}
		scene.QueryIntersectionInfo(ray, &isect);
		isect.ComputeScatteringFunction(arena);
		BSDF* pBSDF = isect.mpBSDF;

		path[bounce].mPdfFwd = ConvertSolidToArea(pdfW, prev, vertex);
		path[bounce].mThroughput = throughput;
		path[bounce].mIsDelta = pBSDF->IsDelta();

		Vec3 wi;
		Vec3 f = pBSDF->Sample(-ray.mDir, &wi, &pdfW, Vec3(rand(), rand(), rand()));
		Vec3 estimation;

		estimation = f * std::abs(Dot(isect.mNormal, wi)) / pdfW;
		throughput = throughput * estimation;

		ray = isect.SpawnRay(wi);

		real pdfWPrev = isect.mpBSDF->Pdf(wi, -1 * ray.mDir);
		prev.mPdfPrev = ConvertSolidToArea(pdfWPrev, vertex, prev);

		bounce++;
		if (bounce >= mMaxDepth + 1) break;
	}

	return bounce;
}

Vec3 PathTracingV2::ConnectToLight(const Scene& scene, StateSequence& rand, const std::vector<PathVertex>& path, int t)
{
	Vec3 L(0, 0, 0);
	const PathVertex& vertex = path[t];
	if (t > 0 && path[t - 1].mIsDelta && vertex.mIsect.mpPrimitive->IsLight()) {
		const Light* pLight = vertex.mIsect.mpPrimitive->GetLight();
		L = vertex.mThroughput * pLight->Emission();
	}
	else {
		L = vertex.mThroughput * SimpleDirectIllumination(scene, vertex.mIsect, rand);
	}
	return L;
}

Vec3 PathTracingV2::SimpleDirectIllumination(const Scene& scene, const Intersection& hitPoint, StateSequence& rand) const
{
	Vec3 L(0, 0, 0);
	if(!hitPoint.mIsDelta){
		Light* pLight;
		real pdfLight;
		real pdfA, pdfW;
		Intersection lightPoint;
		pLight = scene.SampleOneLight(&pdfLight, rand());
		Vec3 Le = pLight->SampleOnePoint(&lightPoint, &pdfA, Vec2(rand(), rand()));
		Vec3 hitToLight = lightPoint.mPos - hitPoint.mPos;
		real dis = hitToLight.Length();
		hitToLight.Normalize();
		real cosTheta0 = hitToLight.Dot(hitPoint.mNormal);
		real cosTheta1 = (-1 * hitToLight).Dot(lightPoint.mNormal);
		pdfW = pdfA * dis * dis / std::abs(cosTheta1);
		Vec3 f = hitPoint.mpBSDF->Evaluate(hitPoint.mOutDir, hitToLight);
		VisibilityTester visibilityTester(hitPoint, lightPoint);
		if (!visibilityTester.Unoccluded(scene) || cosTheta1 < 0) {
			return Vec3(0, 0, 0);
		}
		else {
			L = Le * f * std::abs(cosTheta0) / pdfW / pdfLight;
		}
	}
	return L;
}

Vec3 PathTracingV2::Li(const Scene& scene, StateSequence& rand, MemoryPool& arena, const Ray& r)
{
	Ray ray = r;
	bool deltaBoundEvent = false;
	Vec3 throughput = Vec3(1, 1, 1);
	Vec3 L;
	int i = 1;
	while(true) {
	//for (int i = 0; i < mMaxDepth; ++i) {
		Intersection isect;
		if (!scene.Intersect(ray, &isect)) {
			break;
		}
		scene.QueryIntersectionInfo(ray, &isect);
		isect.ComputeScatteringFunction(arena);
		BSDF* pBSDF = isect.mpBSDF;

		if ((i == 0 || deltaBoundEvent) && isect.mpPrimitive->IsLight()) {
			//GYT_Print("Delta!\n");
			const Light* pLight = isect.mpPrimitive->GetLight();
			L += throughput * pLight->Emission();
		}
		else {
			L += throughput * SimpleDirectIllumination(scene, isect, rand);
		}

		Vec3 wi;
		real pdfW;
		Vec3 f = pBSDF->Sample(-ray.mDir, &wi, &pdfW, Vec3(rand(), rand(), rand()));
		if (f == Vec3() || pdfW == 0) break;
		Vec3 estimation = f * std::abs(Dot(isect.mNormal, wi)) / pdfW;
		deltaBoundEvent = pBSDF->IsDelta();

		real p = std::min((real)1.0, (estimation * throughput).Y() / throughput.Y());
		if (p < 1 && i > 5) {
			if (rand() < p) {
				throughput = throughput / p;
			}
			else {
				break;
			}
		}
		throughput = throughput * estimation;
		ray = isect.SpawnRay(wi);
		i++;
		if (i >= mMaxDepth + 1) break;
	}
	return L;
}

int PathTracingV2::ConstructCameraPath(const Scene& scene, const Camera& camera, StateSequence& rand, MemoryPool& arena, const Ray& cameraRay, std::vector<PathVertex>& cameraPath, int maxDepth)
{
	if (maxDepth == 0) return 0;
	Vec3 throughput(1.f, 1.f, 1.f);
	cameraPath[0].mIsect.mPos = camera.mPos;
	cameraPath[0].mIsect.mOutDir = -1 * cameraRay.mDir;
	cameraPath[0].mIsect.mNormal = camera.mCz;
	cameraPath[0].mIsect.mGeometryNormal = camera.mCz;
	cameraPath[0].mPdfFwd = camera.PdfPos();
	cameraPath[0].mThroughput = throughput;
	real pdfDir = camera.PdfDir(cameraRay);
	const Camera* pCam = &camera;
	cameraPath[0].mpCamera = const_cast<Camera*>(pCam);
	Ray ray = cameraRay;

	real pdfW = pdfDir;
	int i = 1;
	while(true){
		Intersection &isect = cameraPath[i].mIsect;
		if (!scene.Intersect(ray, &isect)) {
			break;
		}
		scene.QueryIntersectionInfo(ray, &isect);
		isect.ComputeScatteringFunction(arena);
		BSDF* pBSDF = isect.mpBSDF;

		cameraPath[i].mPdfFwd = pdfW;
		cameraPath[i].mThroughput = throughput;
		cameraPath[i].mIsDelta = pBSDF->IsDelta();

		//if ((i == 0 || deltaBoundEvent) && isect.mpPrimitive->IsLight()) {
		//	std::cout << i << std::endl;
		//	const Light* pLight = isect.mpPrimitive->GetLight();
		//	L += throughput * pLight->Emission();
		//}
		//else {
		//	L += throughput * SimpleDirectIllumination(scene, isect, rand);
		//}

		Vec3 wi;
		Vec3 f = pBSDF->Sample(-ray.mDir, &wi, &pdfW, Vec3(rand(), rand(), rand()));
		Vec3 estimation;

		estimation = f * std::abs(Dot(isect.mNormal, wi)) / pdfW;
		throughput = throughput * estimation;

		ray = isect.SpawnRay(wi);
		i++;
		if (i >= mMaxDepth + 1) break;
	}

	return i;
}

GYT_NAMESPACE_END