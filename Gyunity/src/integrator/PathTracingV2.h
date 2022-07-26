#pragma once

#include "Integrator.h"

GYT_NAMESPACE_BEGIN

class StateSequence;
class PathTracingV2 : public Integrator {
private:
	int mSpp , mMaxDepth;
	std::shared_ptr<Sampler> mpSampler;
	std::shared_ptr<SamplerEnum> mpSamplerEnum;
public:
	PathTracingV2(const std::shared_ptr<Sampler>& pSampler, const std::shared_ptr<SamplerEnum>& pSamplerEnum, int spp, int maxDepth)
		: mSpp(spp), mMaxDepth(maxDepth), mpSampler(pSampler), mpSamplerEnum(pSamplerEnum) {}
	void Render(const Scene& scene, const Camera& camera);

private:
	Vec3 Li(const Scene& scene, StateSequence& rand, MemoryPool& arena, const Ray& r);

	int ConstructCameraPath(const Scene& scene, const Camera& camera, StateSequence& rand, MemoryPool& arena, const Ray& r, std::vector<PathVertex>& cameraPath, int maxDepth);

	int GenerateCameraPath(
		const Scene				&scene,
		const Camera			&camera,
		StateSequence			&rand,
		MemoryPool				&arena,
		std::vector<PathVertex>	&cameraPath,
		const Ray				&cameraRay,
		int						maxdepth);

	int Trace(
		const Scene				&scene,
		MemoryPool				&arena,
		StateSequence			&rand,
		const Ray				&r,
		int						depth,
		Vec3					throughput,
		real					pdfFwd,
		std::vector<PathVertex>	&lightPath,
		int						maxDepth,
		TransportMode			mode);

	Vec3 ConnectToLight(const Scene& scene, StateSequence& rand, const std::vector<PathVertex>& path, int t);

	Vec3 SimpleDirectIllumination(const Scene& scene, const Intersection& hitPoint, StateSequence& rand) const;
};

GYT_NAMESPACE_END