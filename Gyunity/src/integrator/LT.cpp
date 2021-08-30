
#include "LT.h"
#include "PathVertex.h"

GYT_NAMESPACE_BEGIN


void LightTracing::Render(const Scene& scene, const Camera& camera)
{
	RandomStateSequence rand(mpSampler, 123);
}

void LightTracing::GenerateLightPath(const Scene &scene, const Camera &camera, StateSequence& rand, std::vector<PathVertex>& lightPath, int maxDepth)
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
		light->SampleOnLight(&isect, &dir, &pdfA, &pdfDir, rand(), rand());
	}
	else if (light->IsDeltaLight()) {
		Le = light->Emission();
		light->SampleLight(&isect, &dir, &pdfA, &pdfDir, rand(), rand());
	}
	real cosTheta = isect.mNormal.Dot(dir);
	lightPath[0].isect.mPos = isect.mPos;
	lightPath[0].isect.mNormal = isect.mNormal;
	lightPath[0].Throughput = Le; // / Pdfdir / Pdfpos * CosTheta;
	lightPath[0].PdfFwd = pdfA;
	lightPath[0].isect.mIsDelta = light->IsDeltaLight();
	Vec3 Throughput = lightPath[0].Throughput * cosTheta / lightPdf / pdfA / pdfDir;
	Ray ray(isect.mPos + dir * RayEps, dir);

	Trace(ray, Throughput, pdfDir, rand, lightPath, 1, maxDepth, scene, camera);

	return;
}

int LightTracing::Trace(const Ray& ray, Vec3 throughput, real pdfFwd, StateSequence& rand, std::vector<PathVertex>& lightPath, int depth, int maxDepth, const Scene& scene, const Camera& camera)
{
	real pdfdir = pdfFwd;
	Ray r = ray;
	int bound = depth;
	real pdfW = pdfFwd;

	while (1) {
		Intersection& isect = lightPath[bound].isect;
		real t;
		int id;
		if (!scene.Intersect(r, &isect)) break;

		lightPath[bound].Throughput = throughput;
		lightPath[bound].PdfFwd = pdfW;

		++bound;
		if (bound >= maxDepth + 1) break;

		Vec3 wo;
		Vec3 f = isect.mpBSDF->Sample(-1 * r.mDir, &wo, &pdfW, rand());
		Normalize(wo);
		throughput = throughput * f * (wo.Dot(isect.mNormal)) / pdfW;

		r.mOrig = isect.mPos + wo * RayEps;
		r.mDir = wo;
	}
}



Vec3 LightTracing::WorldToScreen(const Camera & camera, const Vec3& vertex, bool* isInScreen) const
{
	*isInScreen = false;
	Transform WorldToRaster = camera.CameraToRaster * camera.WorldToCamera;
	Vec3 rasterPos = WorldToRaster(vertex);
	if (rasterPos.x >= 0 && rasterPos.x <= camera.GetFilm()->mWidth &&
		rasterPos.y >= 0 && rasterPos.y <= camera.GetFilm()->mHeight) {
		*isInScreen = true;
	}
	return rasterPos;
}


Vec3 LightTracing::ConnectToCamera(const PathVertex& vertex, int s, const Scene& scene, const Camera& camera, StateSequence& rand, Vec3* pRaster, bool* inScreen) {
	if (vertex.isect.mIsDelta) return Vec3(0.0, 0.0, 0.0);

	*inScreen = true;
	Vec3 hitPointToCam = camera.mPos - vertex.isect.mPos;
	double ray_tmax = hitPointToCam.Length();
	hitPointToCam.Norm();
	if (camera.mCz.Dot(-1 * hitPointToCam) < 0) {
		*inScreen = false;
		return Vec3(0.0, 0.0, 0.0);
	}

	bool isInScreen;
	*pRaster = WorldToScreen(camera, vertex.isect.mPos, &isInScreen);
	if (!isInScreen) {
		*inScreen = false;
		return Vec3(0.0, 0.0, 0.0);
	}

	Ray ray(vertex.isect.mPos, hitPointToCam, 0.f, ray_tmax);
	Intersection isc;
	double t; int id;

	if (scene.Intersect(ray)) {
		* inScreen = false;
		return Vec3(0.0, 0.0, 0.0);
	}
	if (s == 0) {
		return vertex.Throughput; //see the light source directly
	}
	real pdfW;
	Vec3 wi;
	Vec3 We = camera.Sample_Wi(vertex.isect, &pdfW, &wi, rand());
	Vec3 f = vertex.isect.mpBSDF->Evaluate(vertex.isect.mOutDir, wi);
	Vec3 L = We * vertex.Throughput * f * wi.Dot(vertex.isect.mNormal) / pdfW;
	return L;
	}


GYT_NAMESPACE_END
