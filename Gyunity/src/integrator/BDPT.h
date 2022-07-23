#pragma once

#include "Integrator.h"
#include "ForwardDecl.h"
#include "PathVertex.h"

GYT_NAMESPACE_BEGIN

class BDPT : public Integrator {
public:
	BDPT() {}
	BDPT(std::shared_ptr<Sampler> pSampler, std::shared_ptr<SamplerEnum> pSamplerEnum, int maxDepth, int spp, bool visualizeStrategies, bool visualizeWeight) :
		mpSampler(pSampler), mpSamplerEnum(pSamplerEnum), mMaxDepth(maxDepth), mSpp(spp), mVisualizeStrategies(visualizeStrategies), mVisualizeWeight(visualizeWeight) {}
	void Render(const Scene& scene, const Camera& camera);
private:
	bool mVisualizeStrategies;
	bool mVisualizeWeight;
	std::shared_ptr<Sampler> mpSampler;
	std::shared_ptr<SamplerEnum> mpSamplerEnum;
	int mSpp, mMaxDepth;

// Private method
	Vec3 ConnectToCamera(const Scene& scene, const Camera& camera, StateSequence& rand,const PathVertex& vertex, int s, Vec3* pRaster, bool* inScreen) const;
	Vec3 ConnectLightToCamera(const Scene& scene, const Camera& camera, const PathVertex& lightVertex, const PathVertex& cameraVertex, Vec2* pRaster);
	Vec3 WorldToScreen(const Camera& camera, const Vec3& point, bool* inScreen) const;
	Vec3 ConnectToLight(const Scene& scene, StateSequence& rand, const std::vector<PathVertex>& path, int t);
};

GYT_NAMESPACE_END