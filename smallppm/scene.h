#pragma once

#include "camera.h"
#include "shape.h"
#include "light.h"
#include "sampling.h"
#include "accelerator.h"
#include "debug_utils.h"
#include "mesh.h"


NAMESPACE_BEGIN

class Scene {
public:
	void Initialize() {
		for (auto light : lights) {
			light->Initialize(*this);
		}
		lightPowerDistribution = ComputeLightPowerDistribution();
		accelerator->SetPrimitives(primitives);
		shapeNum = 0;
	}

	void SetAccelerator(const std::shared_ptr<Accelerator>& pAccelerator) {
		accelerator = pAccelerator;
	}

	void AddPrimitive(std::shared_ptr<Shape> shape, std::shared_ptr<Material> material, MediumInterface mi = MediumInterface()) {
		shape->shapeId = shapeNum;
		std::shared_ptr<Primitive> primitive =
			std::make_shared<GeometryPrimitive>(shape, material, nullptr, mi);
		primitive->primId = primitiveNum;
		primitives.push_back(primitive);
		++shapeNum;
		++primitiveNum;
	}

	void AddPrimitive(std::shared_ptr<Primitive> primitive) {
		Shape* shape = primitive->GetShape();
		if (shape) {
			shape->shapeId = shapeNum;
			++shapeNum;
		}
		primitive->primId = primitiveNum;
		primitives.push_back(primitive);
		++primitiveNum;
	}

	void AddLight(std::shared_ptr<Light> light) {
		lights.push_back(light);
		if (light->IsEnvironmentLight()) {
			envLight = light;
		}
	}


	void AddMesh(Mesh &mesh, const Transform &transform) {
		for (int i = 0; i < mesh.untransformedTriangles.size(); ++i) {
			Triangle* triangle = new Triangle();
			*triangle = mesh.untransformedTriangles[i];
			triangle->SetTransform(transform);
			std::shared_ptr<Shape> shape = std::shared_ptr<Shape>(triangle);
			AddPrimitive(shape, mesh.material, mesh.mediumInterface);
		}
	}

	const Light* GetEnvironmentLight() const {
		return envLight.get();
	}

	bool Intersect(const Ray& r, Intersection* isect) const {
		return accelerator->Intersect(r, isect);
	}

	bool Intersect(const Ray& r) const {
		return accelerator->Intersect(r);
	}

	bool IntersectTr(Ray& ray, StateSequence& rand, Intersection* isect, Vec3* Tr) const;

	const std::vector<std::shared_ptr<Light>>& GetLights() const {
		return lights;
	}

	Light* SampleOneLight(real *lightPdf, real u) const {
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
		return lights[lightNum].get();
	}

	const std::vector<std::shared_ptr<Primitive>>& GetPrimitives() const {
		return primitives;
	}

	void QueryIntersectionInfo(const Ray& ray, Intersection* isect) const;

	void GetBoundingSphere(Vec3* center, real* radius) const;

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
	std::vector<std::shared_ptr<Primitive>> primitives;
	std::vector<std::shared_ptr<Light>> lights;
	std::shared_ptr<Light> envLight = nullptr;
	std::unique_ptr<Distribution1D> lightPowerDistribution;
	std::shared_ptr<Accelerator> accelerator;
	int64 shapeNum;
	int64 primitiveNum;
};

NAMESPACE_END