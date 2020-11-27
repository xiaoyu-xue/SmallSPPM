#include "path.h"
#include "sampler/sampler.h"
#include "system/memory.h"
#include "common/debug_utils.h"
#include "visual/intersection.h"
#include <algorithm>


NAMESPACE_BEGIN

Vec3 PathTracing::Li(const Ray &r, const Scene& scene, StateSequence& rand, MemoryPool& arena) const {
	Ray ray = r;
	bool deltaBoundEvent = false;
	Vec3 throughput = Vec3(1, 1, 1);
	Vec3 L;
	for (int i = 0; i < maxDepth; ++i) {
		Intersection isect;
		if (!scene.Intersect(ray, &isect)) {
			if ((i == 0 || deltaBoundEvent) && scene.GetEnvironmentLight()) {
				return throughput * scene.GetEnvironmentLight()->Emission(ray);
			}
			break;
		}
		scene.QueryIntersectionInfo(ray, &isect);
		isect.ComputeScatteringFunction(arena);
		BSDF* bsdf = isect.bsdf;

		DEBUG_PIXEL_IF(ThreadIndex()) {
			std::cout << "Depth: " << i << " ************************************ \n";
			std::cout << "wo: " << -ray.d << std::endl;
		}

		if ((i == 0 || deltaBoundEvent) && isect.primitive->IsLight()) {
			const Light* emissionShape = isect.primitive->GetLight();
			L += throughput * emissionShape->Emission(isect, isect.wo);
		}
		else {
			L += throughput * DirectIllumination(scene, isect, rand(), Vec2(rand(), rand()), Vec3(rand(), rand(), rand()), rand);

		}

		DEBUG_PIXEL_IF(ThreadIndex()) {
			std::cout << "Path Tracing: Sample BSDF " << std::endl;
		}
		Vec3 wi;
		real pdf;
		Vec3 f = bsdf->Sample_f(-ray.d, &wi, &pdf, Vec3(rand(), rand(), rand()));
		if (f == Vec3() || pdf == 0) break;
		Vec3 estimation = f * std::abs(Dot(isect.n, wi)) / pdf;
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

NAMESPACE_END