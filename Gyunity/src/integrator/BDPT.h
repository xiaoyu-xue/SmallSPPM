#pragma once

#include "Integrator.h"
#include "ForwardDecl.h"
#include "PathVertex.h"

GYT_NAMESPACE_BEGIN

class BDPT : public Integrator {
public:
	BDPT() {}
	BDPT(std::shared_ptr<Sampler> pSampler, int maxDepth, int spp, bool visualizeStrategies, bool visualizeWeight) :
		mpSampler(pSampler), mMaxDepth(maxDepth), mSpp(spp), mVisualizeStrategies(visualizeStrategies), mVisualizeWeight(visualizeWeight) {}
	void Render(const Scene& scene, const Camera& camera);
private:
	bool mVisualizeStrategies;
	bool mVisualizeWeight;
	std::shared_ptr<Sampler> mpSampler;
	int mSpp, mMaxDepth;

// Debug
	Vec3 ConnectToCamera(const Scene& scene, const Camera &camera, StateSequence& rand, const std::vector<PathVertex> &path, int s, Vec2 *pRaster);
	Vec3 ConnectLightToCamera(const Scene& scene, const Camera& camera, const PathVertex& lightVertex, const PathVertex& cameraVertex, Vec2* pRaster);
};

GYT_NAMESPACE_END