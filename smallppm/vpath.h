#pragma once
#include "tiled_integrator.h"
class VolPathTracing : public TiledIntegrator {
public:
	VolPathTracing(int samplePerPixel, int maxDepth, const std::shared_ptr<Sampler>& pSampler, const std::shared_ptr<SamplerEnum>& pSamplerEnum) :
		TiledIntegrator(samplePerPixel, pSampler, pSamplerEnum), maxDepth(maxDepth)
	{

	}

	Vec3 Li(const Ray& ray, const Scene& scene, StateSequence& rand, MemoryArena& arena) const override;

private:
	int maxDepth;
};