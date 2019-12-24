#pragma once
#include "texture.h"

NAMESPACE_BEGIN

class ConstantTexture : public Texture {
public:
	ConstantTexture(Vec3 value) : value(value) {}

	Vec3 Sample(const Vec3 &coord) const override {
		return value;
	}

	Vec3 Sample(const Vec2 &coord) const override {
		return value;
	}

	Vec3 Sample(const Intersection &coord) const override {
		return value;
	}

private:
	Vec3 value;
};

NAMESPACE_END