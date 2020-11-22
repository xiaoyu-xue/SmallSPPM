#pragma once

#include "material.h"
#include "visual/microfacet.h"
#include "texture/texture.h"
#include "bsdf/bsdf.h"

NAMESPACE_BEGIN

class RoughDielectricMaterial : public Material {
public:
    RoughDielectricMaterial(const std::shared_ptr<Texture<Vec3>>& Kr,
        const std::shared_ptr<Texture<Vec3>>& Kt,
        const std::shared_ptr<Texture<real>>& xRoughness,
        const std::shared_ptr<Texture<real>>& yRoughness,
        const std::shared_ptr<Texture<real>>& index) :
        Kr(Kr), Kt(Kt), xRoughness(xRoughness), yRoughness(yRoughness), index(index) {

    }

	void ComputeScatteringFunction(Intersection* isect, MemoryArena& arena, TransportMode mode = TransportMode::Radiance) const {
		real alphax = GGXDistribution::RoughnessToAlpha(xRoughness->Sample(*isect));
		real alphay = GGXDistribution::RoughnessToAlpha(yRoughness->Sample(*isect));
		//real alphax = xRoughness->Sample(*isect);
		//real alphay = yRoughness->Sample(*isect);
		isect->bsdf = ARENA_ALLOC(arena, BSDF)(*isect);

		Vec3 kr = Kr->Sample(*isect);
		Vec3 kt = Kt->Sample(*isect);
		real eta = index->Sample(*isect);
		Fresnel* fresnel = ARENA_ALLOC(arena, FresnelDielectric)(1.f, eta);
		MicrofacetDistribution* distribution =
			ARENA_ALLOC(arena, GGXDistribution)(alphax, alphay, true);

		BxDF* microfacetReflectionBSDF = ARENA_ALLOC(arena, MicrofacetReflectionBSDF)(distribution, fresnel, kr);
		BxDF* microfacetTransmissionBSDF = ARENA_ALLOC(arena, MicrofacetTransmissionBSDF)(distribution, kt, 1.0, eta, mode);
		//isect->bsdf->Add(microfacetReflectionBSDF);
		//isect->bsdf->Add(microfacetTransmissionBSDF);
		BxDF* roughDielectricBSDF = ARENA_ALLOC(arena, RoughDielectricBSDF)(distribution, kr, kt, 1.f, eta, mode);
		isect->bsdf->Add(roughDielectricBSDF);
	}

private:
    std::shared_ptr<Texture<Vec3>> Kr, Kt;
    std::shared_ptr<Texture<real>> xRoughness, yRoughness;
    std::shared_ptr<Texture<real>> index;
};

NAMESPACE_END