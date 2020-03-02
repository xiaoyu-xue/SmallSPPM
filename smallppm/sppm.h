#pragma once

#include "utils.h"
#include "linagl.h"
#include "integrator.h"
#include "sampler_enum.h"
#include "sampler.h"
#include "hashgrid.h"
#include "debug_utils.h"

NAMESPACE_BEGIN

struct HPoint {
	Vec3 importance, pos, nrm, flux, outDir;
	real r2;
	int64 n; // n = N / ALPHA in the paper
	int64 pix;
	bool used;
	int64 m;
};

class Scene;
class SPPM : public Integrator {
public:
	SPPM(int iterations, int nPhotonsPerStage, int maxDepth, real initialRadius, real alpha, bool batchShrink, const std::shared_ptr<Sampler> &pSampler,
		const std::shared_ptr<SamplerEnum> &pSmplerEnum, bool traceGlossyRay = false) :
		nIterations(iterations), nPhotonsPerRenderStage(nPhotonsPerStage), batchShrink(batchShrink),
		maxDepth(maxDepth), initialRadius(initialRadius), alpha(alpha), sampler(pSampler), samplerEnum(pSmplerEnum), traceGlossyRay(traceGlossyRay)
	{

	}

	void GeneratePhoton(const Scene& scene, Ray* pr, Vec3* f, real u, const Vec2& v, const Vec2& w);

	void TraceEyePath(const Scene& scene, StateSequence& rand, const Ray& ray, int64 pixel, MemoryArena& arena);

	void TracePhoton(const Scene& scene, StateSequence& rand, const Ray& ray, Vec3 photonFlux, MemoryArena& arena);

	void GenerateRadiusImage(const Scene& scene, const Camera &camera);

	void Render(const Scene& scene, const Camera &camera) override;

private:
	void Initialize(int w, int h) {
		radius2.resize(w * h);
		photonNums.resize(w * h);
		flux.resize(w * h);
		hitPoints.resize(w * h);
		directillum.resize(w * h);
		c.resize(w * h);
		pixelLocks.resize(w * h);
		for (int i = 0; i < w * h; ++i) {
			radius2[i] = initialRadius * initialRadius;
			photonNums[i] = 0;
			flux[i] = Vec3();
			hitPoints[i].m = 0;
		}
	}

	HashGrid<HPoint*> hashGrid;
	std::shared_ptr<Sampler> sampler;
	std::shared_ptr<SamplerEnum> samplerEnum;
	std::vector<real> radius2;
	//std::vector<int64> photonNums;
	std::vector<real> photonNums;
	std::vector<Vec3> flux;
	std::vector<Vec3> directillum;
	std::vector<HPoint> hitPoints;
	std::vector<Vec3> c;
	std::vector<Spinlock> pixelLocks;
	const real initialRadius;
	const int maxDepth;
	const int nIterations;
	const int64 nPhotonsPerRenderStage;
	const real alpha;
	bool batchShrink;
	bool traceGlossyRay;
};

NAMESPACE_END