#include "vpath.h"
#include "sampler.h"
#include <algorithm>

NAMESPACE_BEGIN

Vec3 VolPathTracing::Li(const Ray& r, const Scene& scene, StateSequence& rand, MemoryArena& arena) const {
	Ray ray = r;
	bool deltaBoundEvent = false;
	Vec3 throughput = Vec3(1, 1, 1);
	Vec3 L;
	for (int i = 0; i < maxDepth; ++i) {
		Intersection isect;
		bool intersect = scene.Intersect(ray, &isect);
		MediumIntersection mi;
		Intersection eqLightPoint;
		if (useEquiAngularSample) {

		}
		else {
			if (ray.medium) throughput *= ray.medium->Sample(ray, rand, arena, &mi);
		}

		if (throughput == Vec3()) break;

		if (mi.IsValid()) {
			if(intersect) scene.QueryIntersectionInfo(ray, &isect);
			L += throughput * DirectIllumination(scene, mi, rand(), Vec2(rand(), rand()), Vec3(rand(), rand(), rand()), rand, true);
			Vec3 wi;
			mi.phase->Sample_p(-ray.d, &wi, Vec2(rand(), rand()));
			ray = mi.SpawnRay(wi);
			deltaBoundEvent = false;
		}
		else {
			if (!intersect) {
				if ((i == 0 || deltaBoundEvent) && scene.GetEnvironmentLight()) {
					return throughput * scene.GetEnvironmentLight()->Emission(ray);
				}
				break;
			}
			scene.QueryIntersectionInfo(ray, &isect);
			isect.ComputeScatteringFunction(arena);
			BSDF* bsdf = isect.bsdf;
			if (!bsdf) {
				ray = isect.SpawnRay(ray.d);
				i--;
				continue;
			}

			if ((i == 0 || deltaBoundEvent) && isect.primitive->IsLight()) {
				std::shared_ptr<Light> emissionShape = isect.primitive->GetLight();
				L += throughput * emissionShape->Emission(isect, isect.wo);
			}
			else {
				L += throughput * DirectIllumination(scene, isect, rand(), Vec2(rand(), rand()), Vec3(rand(), rand(), rand()), rand, true);

			}

			Vec3 wi;
			real pdf;
			Vec3 f = bsdf->Sample_f(-ray.d, &wi, &pdf, Vec3(rand(), rand(), rand()));
			if (f == Vec3() || pdf == 0) break;
			Vec3 estimation = f * std::abs(Dot(isect.n, wi)) / pdf;
			deltaBoundEvent = bsdf->IsDelta();

			throughput = throughput * estimation;
			ray = isect.SpawnRay(wi);
		}
		real p = std::min((real)1.0, (throughput).Y() / throughput.Y());
		if (p < 1 && i > 5) {
			if (rand() < p) {
				throughput = throughput / p;
			}
			else {
				break;
			}
		}
	}
	return L;
}

NAMESPACE_END