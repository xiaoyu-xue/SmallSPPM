#pragma once

#include "camera.h"
#include "shape.h"
#include "light.h"
#include "sampling.h"
#include "accelerator.h"
#include "debug_utils.h"


NAMESPACE_BEGIN

class Scene {
public:

	void SetCamera(const std::shared_ptr<Camera> &pCamera) {
		camera = pCamera;
		shapeNum = 0;
	}

	void SetAccelerator(const std::shared_ptr<Accelerator>& pAccelerator) {
		accelerator = pAccelerator;
	}

	void AddShape(std::shared_ptr<Shape> shape) {
		shapes.push_back(shape);
		shapes[shapeNum]->shapeId = shapeNum;
		++shapeNum;
	}

	void AddPrimitive(std::shared_ptr<Shape> shape, std::shared_ptr<Material> material) {
		shape->shapeId = shapeNum;
		std::shared_ptr<Primitive> primitive =
			std::make_shared<GeometryPrimitive>(shape, material);
		primitives.push_back(primitive);
		++shapeNum;
	}

	void AddPrimitive(std::shared_ptr<Primitive> primitive) {
		std::shared_ptr<Shape> shape = primitive->GetShape();
		if (shape) {
			shape->shapeId = shapeNum;
#ifdef _DEBUG
			primitive->primitiveId = shapeNum;
#endif
			++shapeNum;
		}
		primitives.push_back(primitive);
	}

	void AddLight(std::shared_ptr<Light> light) {
		lights.push_back(light);
		if (light->IsAreaLight()) {
			AddShape(light->GetShapePtr());
		}
	}

	void AddLight(std::shared_ptr<Primitive> lightPrimitive) {
		std::shared_ptr<Light> light = lightPrimitive->GetLight();
		lights.push_back(light);
		if (light->IsAreaLight()) {
			AddPrimitive(lightPrimitive);
		}
	}

	void Initialize() {
		lightPowerDistribution = ComputeLightPowerDistribution();
		accelerator->SetPrimitives(primitives);
	}

	//bool Intersect(const Ray &r, real *t, Intersection *isect, std::shared_ptr<Shape> &hitObj) const {
	//	int n = (int)shapes.size();
	//	*t = Inf;
	//	for (int i = 0; i < n; ++i) {
	//		Intersection intersection;
	//		real ti;
	//		if (shapes[i]->Intersect(r, &intersection, &ti)) {
	//			if (ti < *t) {
	//				*t = ti;
	//				*isect = intersection;
	//				hitObj = shapes[i];
	//			}
	//		}
	//	}
	//	return  (*t > r.tMin && *t < r.tMax);
	//}

	//bool Intersect(const Ray &r) const {
	//	for (auto shape : shapes) {
	//		if (shape->Intersect(r)) {
	//			return true;
	//		}
	//	}
	//	return false;
	//}

	bool Intersect(const Ray& r, real* t, Intersection* isect) const {
		return accelerator->Intersect(r, t, isect);
	}

	bool Intersect(const Ray& r) const {
		return accelerator->Intersect(r);
	}

	const std::vector<std::shared_ptr<Light>>& GetLights() const {
		return lights;
	}

	const std::vector<std::shared_ptr<Shape>>& GetShapes() const {
		return shapes;
	}

	std::shared_ptr<Camera> GetCamera() const {
		return camera;
	}

	std::shared_ptr<Light> SampleOneLight(real *lightPdf, real u) const {
		int nLights = (int)(lights.size());
		int lightNum;
		if (lightPowerDistribution != nullptr) {
			lightNum = lightPowerDistribution->SampleDiscrete(u, lightPdf);
			if (*lightPdf == 0) return nullptr;
		}
		else {
			lightNum = std::min((int)(u * nLights), nLights - 1);
			*lightPdf = 1.f / nLights;
		}
		return lights[lightNum];
	}

private:

	std::unique_ptr<Distribution1D> ComputeLightPowerDistribution() {
		if (lights.size() == 0) return nullptr;
		std::vector<real> lightPower;
		for (const auto &light : lights)
			lightPower.push_back(light->Power().Y());
		return std::unique_ptr<Distribution1D>(
			new Distribution1D(&lightPower[0], (int)lightPower.size()));
	}

private:
	std::vector<std::shared_ptr<Shape>> shapes;
	std::vector<std::shared_ptr<Primitive>> primitives;
	std::vector<std::shared_ptr<Light>> lights;
	std::shared_ptr<Camera> camera;
	std::unique_ptr<Distribution1D> lightPowerDistribution;
	std::shared_ptr<Accelerator> accelerator;
	int shapeNum;
	int primitiveNum;
};

NAMESPACE_END