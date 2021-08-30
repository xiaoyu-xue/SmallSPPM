#pragma once

#include "ForwardDecl.h"
#include <vector>
#include "Integrator.h"
#include "sampler/Sampler.h"

GYT_NAMESPACE_BEGIN

class PathVertex;
class StateSequence;

class LightTracing : public Integrator
{
public:
	LightTracing(const std::shared_ptr<Sampler> pSampler,  Camera* pCamera, int maxDepth, int spp) :
		spp(spp), mpSampler(pSampler), mpCamera(pCamera), maxDepth(maxDepth), mLightPath(maxDepth + 1){}
	void Render(const Scene& scene, const Camera& camera) override;
private:

	void GenerateLightPath(const Scene& scene, const Camera& camera, StateSequence& rand, std::vector<PathVertex>& lightPath, int maxDepth);
	int Trace(const Ray& ray, Vec3 throughput, real pdfFwd, StateSequence& rand, std::vector<PathVertex>& lightPath, int depth, int maxDepth, const Scene& scene, const Camera& camera);
	Vec3 ConnectToCamera(const PathVertex& vertex, int s, const Scene& scene, const Camera& camera, StateSequence& rand, Vec3* pRaster, bool* inScreen);
	Vec3 WorldToScreen(const Camera &camerea, const Vec3& vertex, bool* isInScreen) const;

	int spp, maxDepth;
	Camera* mpCamera;
	std::shared_ptr<Sampler> mpSampler;
	std::vector<PathVertex> mLightPath;
};


GYT_NAMESPACE_END