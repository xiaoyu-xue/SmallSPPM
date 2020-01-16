#pragma once
#include "def.h"
#include "material.h"
#include "texture.h"
#include "bsdf.h"

NAMESPACE_BEGIN

class MirrorMaterial : public Material {
public:
	MirrorMaterial(const std::shared_ptr<Texture<Vec3>>& kr) : kr(kr) {

	}
	//void ComputeScatteringFunction(Intersection* isect,
	//	TransportMode mode = TransportMode::Radiance) const {

	//	isect->bsdf = std::shared_ptr<BSDF>(new SpecularBSDF(*isect, kr->Sample(*isect)));
	//}
	//void ComputeScatteringFunction(Intersection* isect,
	//	TransportMode mode = TransportMode::Radiance) const {

	//	isect->bsdf = std::shared_ptr<BSDF>(new BSDF(*isect));
	//	std::shared_ptr<BxDF> specularBSDF(new SpecularBSDF(kr->Sample(*isect)));
	//	isect->bsdf->Add(specularBSDF);
	//}

	void ComputeScatteringFunction(Intersection* isect, MemoryArena &arena,
		TransportMode mode = TransportMode::Radiance) const {

		isect->bsdf = ARENA_ALLOC(arena, BSDF)(*isect);

		isect->bsdf->Add(ARENA_ALLOC(arena, SpecularBSDF)(kr->Sample(*isect)));
	}
private:
	std::shared_ptr<Texture<Vec3>> kr;
};

NAMESPACE_END