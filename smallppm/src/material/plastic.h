#pragma once

#include "material.h"
#include "visual/microfacet.h"
#include "texture/texture.h"
#include "bsdf/bsdf.h"

NAMESPACE_BEGIN

class PlasticMaterial : public Material {
public:
	PlasticMaterial(const std::shared_ptr<Texture<Vec3>>& Kd, const std::shared_ptr<Texture<Vec3>>& Ks, const std::shared_ptr<Texture<real>>& roughness):
		kd(Kd), ks(Ks), roughness(roughness){}

	void ComputeScatteringFunction(Intersection* isect, MemoryArena& arena, TransportMode mode = TransportMode::Radiance) const {
		real alpha = GGXDistribution::RoughnessToAlpha(roughness->Sample(*isect));
		isect->bsdf = ARENA_ALLOC(arena, BSDF)(*isect);

		Vec3 Kd = kd->Sample(*isect);
		isect->bsdf->Add(ARENA_ALLOC(arena, DiffuseBSDF)(Kd));
		Fresnel* fresnel = ARENA_ALLOC(arena, FresnelDielectric)(1.5f, 1.f);
		MicrofacetDistribution* distribution =
			ARENA_ALLOC(arena, GGXDistribution)(alpha, alpha, true);
		Vec3 Ks = ks->Sample(*isect);
		BxDF* microfacetReflectionBSDF = ARENA_ALLOC(arena, MicrofacetReflectionBSDF)(distribution, fresnel, Ks);
		isect->bsdf->Add(microfacetReflectionBSDF);
	}
private:
	std::shared_ptr<Texture<Vec3>> kd, ks;
	std::shared_ptr<Texture<real>> roughness;
};

NAMESPACE_END