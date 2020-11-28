#pragma once
#include "Intersection.h"
#include "shape/Shape.h"
#include "material/Material.h"
#include "light/Light.h"
#include "math/AABB.h"
#include "Medium.h"

GY_NAMESPACE_BEGIN

//class Primitive {
//public:
//
//	virtual bool Intersect(const Ray& r, Intersection* isect, real* t) const = 0;
//	virtual bool Intersect(const Ray& r) const = 0;
//	virtual void ComputeScatteringFunction(Intersection* isect, MemoryArena& arena,
//		TransportMode mode = TransportMode::Radiance) const = 0;
//	virtual std::shared_ptr<Shape> GetShape() const = 0;
//	virtual std::shared_ptr<Material> GetMaterial() const = 0;
//	virtual std::shared_ptr<Light> GetLight() const = 0;
//	virtual bool IsLight() const = 0;
//	virtual AABB WorldBound() const = 0;
//	virtual ~Primitive(){}
//
//#ifdef _DEBUG
//	int primitiveId;
//#endif
//};


//class GeometryPrimitive : public Primitive {
//public:
//	GeometryPrimitive(const std::shared_ptr<Shape>& shape, 
//		const std::shared_ptr<Material>& material = nullptr,
//		const std::shared_ptr<Light>& light = nullptr) :
//		shape(shape), material(material), light(light) {
//
//	}
//	bool Intersect(const Ray& r, Intersection* isect, real* t) const override;
//	bool Intersect(const Ray& r) const override;
//	void ComputeScatteringFunction(Intersection* isect, MemoryArena &arena, TransportMode mode = TransportMode::Radiance) const override;
//
//	std::shared_ptr<Shape> GetShape() const override {
//		return shape;
//	}
//
//	std::shared_ptr<Material> GetMaterial() const override {
//		return material;
//	}
//
//	std::shared_ptr<Light> GetLight() const override {
//		return light;
//	}
//
//	bool IsLight() const override {
//		return light != nullptr;
//	}
//
//	AABB WorldBound() const override {
//		return shape->WorldBould();
//	}
//public:
//	std::shared_ptr<Shape> shape;
//	std::shared_ptr<Material> material;
//	std::shared_ptr<Light> light;
//};


//class TransformedPrimitive : public Primitive {
//public:
//
//private:
//	std::shared_ptr<Primitive> primitive;
//};

class  Primitive {
public:
	Primitive(const std::shared_ptr<Shape>& shape,
		const std::shared_ptr<Material>& material = nullptr,
		const std::shared_ptr<Light>& light = nullptr,
		MediumInterface mediumInterface = MediumInterface()) 
		:mpShape(shape), mpMaterial(material), mpLight(light), mMediumInterface(mediumInterface) 
	{

	}
	bool Intersect(const Ray& r, Intersection* isect) const;

	bool Intersect(const Ray& r) const;

	void ComputeScatteringFunction(Intersection* isect, MemoryPool& arena, TransportMode mode = TransportMode::Radiance) const;

	Shape* GetShape() const {
		return mpShape.get();
	}

	const Material* GetMaterial() const {
		return mpMaterial.get();
	}

	const Light* GetLight() const {
		return mpLight.get();
	}

	bool IsLight() const {
		return mpLight != nullptr;
	}

	AABB WorldBound() const {
		return mpShape->WorldBould();
	}

	void QueryIntersectionInfo(const Ray& ray, Intersection* isect) const;

public:
	int64 mPrimId;
	std::shared_ptr<Shape> mpShape;
	std::shared_ptr<Material> mpMaterial;
	std::shared_ptr<Light> mpLight;
	MediumInterface mMediumInterface;

};

using GeometryPrimitive = Primitive;

GY_NAMESPACE_END