#pragma once

#include "ForwardDecl.h"
#include <vector>
#include "Integrator.h"
#include "sampler/Sampler.h"
#include "PathVertex.h"

GYT_NAMESPACE_BEGIN


class StateSequence;

class LightTracing : public Integrator
{
public:
	LightTracing(const std::shared_ptr<Sampler> pSampler, int maxDepth, int spp) :
		spp(spp), mpSampler(pSampler), maxDepth(maxDepth) {}
	
	void Render(const Scene& scene, const Camera& camera) override;

private:

	int GenerateLightPath(const Scene& scene, const Camera& camera, StateSequence& rand, MemoryPool& arena, std::vector<PathVertex> &lightPath);
	int Trace(const Ray& ray, Vec3 throughput, real pdfFwd, StateSequence& rand, int depth, const Scene& scene, const Camera& camera, MemoryPool& arena, std::vector<PathVertex>& lightPath);
	Vec3 ConnectToCamera(const PathVertex& vertex, int s, const Scene& scene, const Camera& camera, StateSequence& rand, Vec3* pRaster, bool* inScreen);
	Vec3 WorldToScreen(const Camera &camerea, const Vec3& vertex, bool* isInScreen) const;

	int spp, maxDepth;
	std::shared_ptr<Sampler> mpSampler;
	//std::vector<PathVertex> mLightPath;

};


GYT_NAMESPACE_END