#pragma once
#include "Texture.h"

GY_NAMESPACE_BEGIN

template<typename T>
class ConstantTexture : public Texture<T> {
public:
	ConstantTexture(T value) : value(value) {}

	T Sample(const Vec3 &coord) const override {
		return value;
	}

	T Sample(const Vec2 &coord) const override {
		return value;
	}

	T Sample(const Intersection &coord) const override {
		return value;
	}

private:
	T value;
};

GY_NAMESPACE_END