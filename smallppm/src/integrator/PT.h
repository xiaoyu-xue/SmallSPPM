#pragma once
#include "TiledIntegrator.h"

class StateSequence;
NAMESPACE_BEGIN


class PathTracing : public TiledIntegrator {
public:
	PathTracing(int samplePerPixel, int maxDepth, const std::shared_ptr<Sampler>& pSampler, const std::shared_ptr<SamplerEnum>& pSamplerEnum):
	TiledIntegrator(samplePerPixel, pSampler, pSamplerEnum), maxDepth(maxDepth)
	{

	}

	Vec3 Li(const Ray &ray, const Scene& scene, StateSequence& rand, MemoryPool& arena) const override;

private:
	int maxDepth;
};

NAMESPACE_END
