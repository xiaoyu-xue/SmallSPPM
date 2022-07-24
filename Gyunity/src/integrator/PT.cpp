#include "PT.h"
#include "common/Core.h"
#include "sampler/Sampler.h"
#include "system/Memory.h"
#include "common/DebugUtils.h"
#include "core/Intersection.h"

GYT_NAMESPACE_BEGIN

Vec3 PathTracing::Li(const Ray& r, const Scene& scene, StateSequence& rand, MemoryPool& arena) const
{
	Ray ray = r;
	bool deltaBoundEvent = false;
	Vec3 throughput = Vec3(1, 1, 1);
	Vec3 L;
	for (int i = 0; i < mMaxDepth; ++i) {
		Intersection isect;
		if (!scene.Intersect(ray, &isect)) {
			if ((i == 0 || deltaBoundEvent) && scene.GetEnvironmentLight()) {
				return throughput * scene.GetEnvironmentLight()->Emission(ray);
			}
			break;
		}
		scene.QueryIntersectionInfo(ray, &isect);
		isect.ComputeScatteringFunction(arena);
		BSDF* bsdf = isect.mpBSDF;

		if ((i == 0 || deltaBoundEvent) && isect.mpPrimitive->IsLight()) {
			const Light* emissionShape = isect.mpPrimitive->GetLight();
			L += throughput * emissionShape->Emission(isect, isect.mOutDir);
		}
		else {
			L += throughput * DirectIllumination(scene, isect, rand(), Vec2(rand(), rand()), Vec3(rand(), rand(), rand()), rand);
		}
		Vec3 wi;
		real pdf;
		Vec3 f = bsdf->Sample(-ray.mDir, &wi, &pdf, Vec3(rand(), rand(), rand()));
		if (f == Vec3() || pdf == 0) break;
		Vec3 estimation = f * std::abs(Dot(isect.mNormal, wi)) / pdf;
		deltaBoundEvent = bsdf->IsDelta();

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

		//if (!(bsdf->IsDelta() || (bsdf->MatchScatterType(ScatterEventType::BSDF_GLOSSY) && i < maxDepth - 1))) break;
	}
	return L;
}

PathTracing::PathTracing(int samplePerPixel, int maxDepth, const std::shared_ptr<Sampler>& pSampler, const std::shared_ptr<SamplerEnum>& pSamplerEnum)
	:TiledIntegrator(samplePerPixel, pSampler, pSamplerEnum), mMaxDepth(maxDepth)
{

}

GYT_NAMESPACE_END