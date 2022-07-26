#include "SimplePT.h"
#include "common/Core.h"
#include "sampler/Sampler.h"
#include "system/Memory.h"
#include "core/Intersection.h"
#include "core/Visibility.h"

GYT_NAMESPACE_BEGIN

Vec3 SimplePathTracing::Li(const Ray& r, const Scene& scene, StateSequence& rand, MemoryPool& arena) const
{
	Ray ray = r;
	bool deltaBoundEvent = false;
	Vec3 throughput = Vec3(1, 1, 1);
	Vec3 L;
	for (int i = 1; i < mMaxDepth; ++i) {
		Intersection isect;
		if (!scene.Intersect(ray, &isect)) {
			break;
		}
		scene.QueryIntersectionInfo(ray, &isect);
		isect.ComputeScatteringFunction(arena);
		BSDF* pBSDF = isect.mpBSDF;

		if ((i == 0 || deltaBoundEvent) && isect.mpPrimitive->IsLight()) {
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
	}
	return L;
}

Vec3 SimplePathTracing::SimpleDirectIllumination(const Scene& scene, const Intersection& hitPoint, StateSequence& rand) const
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

SimplePathTracing::SimplePathTracing(int samplePerPixel, int maxDepth, const std::shared_ptr<Sampler>& pSampler, const std::shared_ptr<SamplerEnum>& pSamplerEnum)
	:TiledIntegrator(samplePerPixel, pSampler, pSamplerEnum), mMaxDepth(maxDepth)
{

}

GYT_NAMESPACE_END