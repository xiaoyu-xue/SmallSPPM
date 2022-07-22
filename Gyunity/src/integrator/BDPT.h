#pragma once

#include "Integrator.h"
#include "ForwardDecl.h"

GYT_NAMESPACE_BEGIN

class BDPT : public Integrator {
public:
	BDPT() {}
	BDPT(std::shared_ptr<SimpleSampler> pSampler, int maxDepth, int spp, bool visualizeStrategies, bool visualizeWeight) :
		mpSampler(pSampler), mMaxDepth(maxDepth), mSpp(spp), mVisualizeStrategies(visualizeStrategies), mVisualizeWeight(visualizeWeight) {}
	void Render(const Scene& scene, const Camera& camera);
private:
	bool mVisualizeStrategies;
	bool mVisualizeWeight;
	std::shared_ptr<SimpleSampler> mpSampler;
	int mSpp, mMaxDepth;
};

GYT_NAMESPACE_END