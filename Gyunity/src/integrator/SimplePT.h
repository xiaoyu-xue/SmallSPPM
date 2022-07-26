#pragma once

#include "TiledIntegrator.h"


GYT_NAMESPACE_BEGIN

class StateSequence;
class SimplePathTracing : public TiledIntegrator {
private:
	int mMaxDepth;
public:
	SimplePathTracing(int samplePerPixel, int maxDepth, const std::shared_ptr<Sampler>& pSampler, const std::shared_ptr<SamplerEnum>& pSamplerEnum);
	Vec3 Li(const Ray& r, const Scene& scene, StateSequence& rand, MemoryPool& arena) const override;
private:
	Vec3 SimpleDirectIllumination(const Scene& scene, const Intersection& isect, StateSequence& rand) const;
};

GYT_NAMESPACE_END