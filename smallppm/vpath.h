#pragma once
#include "tiled_integrator.h"
class VolPathTracing : public TiledIntegrator {
public:
	VolPathTracing(int samplePerPixel, int maxDepth, const std::shared_ptr<Sampler>& pSampler, 
		const std::shared_ptr<SamplerEnum>& pSamplerEnum, bool equalAngular = true) :
		TiledIntegrator(samplePerPixel, pSampler, pSamplerEnum), 
		maxDepth(maxDepth), useEquiAngularSample(equalAngular)
	{

	}

	Vec3 Li(const Ray& ray, const Scene& scene, StateSequence& rand, MemoryArena& arena) const override;

private:
	Vec3 ConnectToLight(const Scene& scene, StateSequence& rand,
		const Intersection& isect, const Intersection& lightPoint,
		const Vec3 &Le, real lightPdf) const;

	int maxDepth;
	bool useEquiAngularSample;
};