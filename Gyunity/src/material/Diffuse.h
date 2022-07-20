#pragma once
#include "common/Core.h"
#include "Material.h"
#include "texture/Texture.h"
#include "bsdf/BSDF.h"
#include "core/Intersection.h"

GYT_NAMESPACE_BEGIN


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
		isect->mpBSDF = MEMORY_POOL_ALLOC(arena, BSDF)(*isect);
		isect->mpBSDF->Add(MEMORY_POOL_ALLOC(arena, DiffuseBSDF)(kd->Sample(*isect)));
	}
private:
	std::shared_ptr<Texture<Vec3>> kd;
};

GYT_NAMESPACE_END