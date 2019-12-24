#pragma once
#include "intersection.h"
#include "shape.h"
#include "material.h"
#include "light.h"

class Primitive {
public:

	virtual bool Intersect(const Ray& r, Intersection* isect, real* t) const = 0;
	virtual bool Intersect(const Ray& r) const = 0;
	virtual void ComputeScatteringFunction(Intersection* isect,
		TransportMode mode = TransportMode::Radiance) = 0;
	virtual std::shared_ptr<Shape> GetShape() const = 0;
	virtual std::shared_ptr<Material> GetMaterial() const = 0;
	virtual std::shared_ptr<Light> GetLight() const = 0;
	virtual bool IsLight() const = 0;
	virtual ~Primitive(){}
};


class GeometryPrimitive : public Primitive {
public:
	GeometryPrimitive(const std::shared_ptr<Shape>& shape, 
		const std::shared_ptr<Material>& material = nullptr,
		const std::shared_ptr<Light>& light = nullptr) :
		shape(shape), material(material), light(light) {

	}
	bool Intersect(const Ray& r, Intersection* isect, real* t) const override;
	bool Intersect(const Ray& r) const override;
	void ComputeScatteringFunction(Intersection* isect, TransportMode mode = TransportMode::Radiance) override;

	std::shared_ptr<Shape> GetShape() const override {
		return shape;
	}

	std::shared_ptr<Material> GetMaterial() const override {
		return material;
	}

	std::shared_ptr<Light> GetLight() const override {
		return light;
	}

	bool IsLight() const override {
		return light != nullptr;
	}
public:
	std::shared_ptr<Shape> shape;
	std::shared_ptr<Material> material;
	std::shared_ptr<Light> light;
};