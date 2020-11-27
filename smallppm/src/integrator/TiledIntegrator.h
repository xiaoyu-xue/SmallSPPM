#pragma once
#include "Integrator.h"



NAMESPACE_BEGIN

class Sampler;
class SamplerEnum;
class StateSequence;

struct Tile {
	int minX, minY, maxX, maxY;
	Tile() {}
	Tile(int minX, int minY, int maxX, int maxY): minX(minX), minY(minY), maxX(maxX), maxY(maxY) {}
};


class TiledIntegrator : public Integrator {
public:
	TiledIntegrator(){}
	TiledIntegrator(int samplePerPixel, const std::shared_ptr<Sampler> &pSampler, const std::shared_ptr<SamplerEnum> &pSamplerEnum): 
		spp(samplePerPixel), sampler(pSampler), samplerEnum(pSamplerEnum) {
	}
	virtual ~TiledIntegrator(){}
	virtual void Render(const Scene& scene, const Camera& camera) override;
	virtual Vec3 Li(const Ray &ray, const Scene& scene, StateSequence &rand, MemoryPool &arena) const = 0;

protected:
	int spp;
	std::shared_ptr<Sampler> sampler;
	std::shared_ptr<SamplerEnum> samplerEnum;

};

NAMESPACE_END
