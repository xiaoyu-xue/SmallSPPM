#pragma once
#include "def.h"
#include "material.h"
#include "texture.h"
#include "bsdf.h"

NAMESPACE_BEGIN

class GlassMaterial : public Material {
public:
	GlassMaterial(const std::shared_ptr<Texture>& kr, const std::shared_ptr<Texture> &kt, real n1 = 1.f, real n2 = 1.5f)
		: kr(kr), kt(kt), eta1(n1), eta2(n2)
	{

	}
	void ComputeScatteringFunction(Intersection* isect,
		TransportMode mode = TransportMode::Radiance) const {

		isect->bsdf =
			std::shared_ptr<BSDF>(new TransmissionBSDF(*isect, kt->Sample(*isect), mode, 1.0, 1.5));
	}
private:
	std::shared_ptr<Texture> kr, kt;
	real eta1, eta2;
};

NAMESPACE_END