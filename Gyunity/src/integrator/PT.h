#pragma once
#include "TiledIntegrator.h"


GY_NAMESPACE_BEGIN

class StateSequence;
class PathTracing : public TiledIntegrator {
private:
	int mMaxDepth;
public:
	PathTracing(int samplePerPixel, int maxDepth, const std::shared_ptr<Sampler>& pSampler, const std::shared_ptr<SamplerEnum>& pSamplerEnum);
	Vec3 Li(const Ray &r, const Scene& scene, StateSequence& rand, MemoryPool& arena) const override;
};

GY_NAMESPACE_END
