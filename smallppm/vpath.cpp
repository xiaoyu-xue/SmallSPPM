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
		real eqLightProb;
		real eqLightPdf;
		Vec3 eqLightLe;
		if (useEquiAngularSample) {
			std::shared_ptr<Light> light = scene.SampleOneLight(&eqLightProb, rand());
			eqLightLe = light->SampleOnePoint(&eqLightPoint, &eqLightPdf, Vec2(rand(), rand()));
			if (ray.medium) throughput *= ray.medium->EquiAngularSampling(ray, rand, arena, eqLightPoint, &mi);
		}
		else {
			if (ray.medium) throughput *= ray.medium->Sample(ray, rand, arena, &mi);
		}

		if (throughput == Vec3()) break;

		if (mi.IsValid()) {
			if (useEquiAngularSample) {
				//if (intersect) scene.QueryIntersectionInfo(ray, &isect);
				//L += throughput * ConnectToLight(scene, rand, mi, eqLightPoint, eqLightLe, eqLightPdf) / eqLightProb;
				//Vec3 wi;
				//mi.phase->Sample_p(-ray.d, &wi, Vec2(rand(), rand()));
				//ray = mi.SpawnRay(wi);
				//deltaBoundEvent = false;
				Vec3 direct = throughput * DirectIllumination(scene, mi, rand(), Vec2(rand(), rand()), Vec3(rand(), rand(), rand()), rand, true);
				if (direct.x < 0 || direct.y < 0 || direct.z < 0) std::cout << i << " " << direct << " " << throughput << std::endl;
				L += direct;
				Vec3 wi;
				mi.phase->Sample_p(-ray.d, &wi, Vec2(rand(), rand()));
				ray = mi.SpawnRay(wi);
				deltaBoundEvent = false;
			}
			else {
				//if (intersect) scene.QueryIntersectionInfo(ray, &isect);
				L += throughput * DirectIllumination(scene, mi, rand(), Vec2(rand(), rand()), Vec3(rand(), rand(), rand()), rand, true);
				Vec3 wi;
				mi.phase->Sample_p(-ray.d, &wi, Vec2(rand(), rand()));
				ray = mi.SpawnRay(wi);
				deltaBoundEvent = false;
			}

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


Vec3 VolPathTracing::ConnectToLight(const Scene &scene, StateSequence &rand,
	const Intersection& isect, const Intersection& lightPoint, const Vec3 &Le, real lightPdf) const {
	Vec3 f;
	real cosTheta = 1;
	Vec3 dir = isect.hit - lightPoint.hit;
	Vec3 wi = -dir.Norm();
	if (lightPoint.n != Vec3()) {
		cosTheta = Dot(dir.Norm(), lightPoint.n);
		if (cosTheta < 0) return Vec3();
	}
	VisibilityTester visibilityTester(isect, lightPoint);
	const MediumIntersection& mi = (const MediumIntersection&)isect;
	real p = mi.phase->p(mi.wo, wi);
	f = Vec3(p);
	return Le * f * visibilityTester.Tr(scene, rand) * cosTheta / dir.Length2() / lightPdf;
}
NAMESPACE_END