#pragma once

#include "def.h"
#include "material.h"
#include "texture.h"
#include "bsdf.h"
#include "microfacet.h"
NAMESPACE_BEGIN

class RoughnessMaterial : public Material {
public:
	RoughnessMaterial(const std::shared_ptr<Texture<Vec3>>& reflectence, 
		const std::shared_ptr<Texture<real>> &roughnessx, const std::shared_ptr<Texture<real>>& roughnessy,
		const std::shared_ptr<Texture<real>>& eta, const std::shared_ptr<Texture<real>>& k)
		: reflectence(reflectence), roughnessx(roughnessx), roughnessy(roughnessy), eta(eta), k(k) {

	}
	//void ComputeScatteringFunction(Intersection* isect,
	//	TransportMode mode = TransportMode::Radiance) const {
	//	MicrofacetDistribution *distribution = new GGXDistribution(alpha);
	//	isect->bsdf = std::shared_ptr<BSDF>(new MicrofacetReflectionBSDF(*isect, distribution, reflectence->Sample(*isect)));
	//}
	//void ComputeScatteringFunction(Intersection* isect,
	//	TransportMode mode = TransportMode::Radiance) const {
	//	real alphax = GGXDistribution::RoughnessToAlpha(roughnessx->Sample(*isect));
	//	real alphay = GGXDistribution::RoughnessToAlpha(roughnessy->Sample(*isect));
	//	MicrofacetDistribution* distribution = 
	//		new GGXDistribution(alphax, alphay, true);
	//	isect->bsdf = std::shared_ptr<BSDF>(new BSDF(*isect));
	//	Fresnel* fresnel = new FresnelConductor(Vec3(1.0, 1.0, 1.0), eta->Sample(*isect), k->Sample(*isect));
	//	std::shared_ptr<BxDF> microfacetReflectionBSDF(new MicrofacetReflectionBSDF(distribution, fresnel, reflectence->Sample(*isect)));
	//	isect->bsdf->Add(microfacetReflectionBSDF);
	//}

	void ComputeScatteringFunction(Intersection* isect, MemoryArena &arena,
		TransportMode mode = TransportMode::Radiance) const {
		real alphax = GGXDistribution::RoughnessToAlpha(roughnessx->Sample(*isect));
		real alphay = GGXDistribution::RoughnessToAlpha(roughnessy->Sample(*isect));
		MicrofacetDistribution* distribution =
			ARENA_ALLOC(arena, GGXDistribution)(alphax, alphay, true);
		isect->bsdf = ARENA_ALLOC(arena, BSDF)(*isect);
		Fresnel* fresnel = ARENA_ALLOC(arena, FresnelConductor)(Vec3(1.0, 1.0, 1.0), eta->Sample(*isect), k->Sample(*isect));
		BxDF* microfacetReflectionBSDF = ARENA_ALLOC(arena, MicrofacetReflectionBSDF)(distribution, fresnel, reflectence->Sample(*isect));
		isect->bsdf->Add(microfacetReflectionBSDF);
	}
private:
	std::shared_ptr<Texture<Vec3>> reflectence;
	std::shared_ptr<Texture<real>> roughnessx, roughnessy;
	std::shared_ptr<Texture<real>> eta;
	std::shared_ptr<Texture<real>> k;
};

NAMESPACE_END