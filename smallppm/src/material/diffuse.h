#pragma once
#include "def.h"
#include "material.h"
#include "texture/texture.h"
#include "bsdf/bsdf.h"

NAMESPACE_BEGIN


class DiffuseMaterial: public Material{
public:
	DiffuseMaterial(const std::shared_ptr<Texture<Vec3>> &kd) : kd(kd) {

	}
	//void ComputeScatteringFunction(Intersection* isect,
	//	TransportMode mode = TransportMode::Radiance) const {
	//	isect->bsdf = std::shared_ptr<BSDF>(new DiffuseBSDF(*isect, kd->Sample(*isect)));
	//}

	//void ComputeScatteringFunction(Intersection* isect,
	//	TransportMode mode = TransportMode::Radiance) const {
	//	isect->bsdf = std::shared_ptr<BSDF>(new BSDF(*isect));
	//	std::shared_ptr<BxDF> diffuseBSDF(new DiffuseBSDF(kd->Sample(*isect)));
	//	isect->bsdf->Add(diffuseBSDF);
	//}
	void ComputeScatteringFunction(Intersection* isect, MemoryPool& arena,
		TransportMode mode = TransportMode::Radiance) const {
		isect->bsdf = MEMORY_POOL_ALLOC(arena, BSDF)(*isect);
		isect->bsdf->Add(MEMORY_POOL_ALLOC(arena, DiffuseBSDF)(kd->Sample(*isect)));
	}
private:
	std::shared_ptr<Texture<Vec3>> kd;
};

NAMESPACE_END