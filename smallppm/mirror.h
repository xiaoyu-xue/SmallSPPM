#pragma once
#include "def.h"
#include "material.h"
#include "texture.h"
#include "bsdf.h"

NAMESPACE_BEGIN

class MirrorMaterial : public Material {
public:
	MirrorMaterial(const std::shared_ptr<Texture>& kr) : kr(kr) {

	}
	void ComputeScatteringFunction(Intersection* isect,
		TransportMode mode = TransportMode::Radiance) const {

		isect->bsdf = std::shared_ptr<BSDF>(new SpecularBSDF(*isect, kr->Sample(*isect)));
	}
private:
	std::shared_ptr<Texture> kr;
};

NAMESPACE_END