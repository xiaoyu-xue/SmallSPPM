#pragma once

#include "light.h"
#include "environment_map.h"

NAMESPACE_BEGIN

class EnvironmentLight : public Light {
public:
	EnvironmentLight(const std::string& filename, real rotate, real scale = 1.f) {
		envMap.reset(new EnvironmentMap(filename, rotate, scale));
	}

	void Initialize(const Scene& scene);

	Vec3 Emission() const override {
		return Vec3();
	}
	Vec3 Emission(const Ray& ray) const override;
	Vec3 SampleLight(Intersection* isect, Vec3* dir, real* pdfPos, real* pdfDir, const Vec2& u, const Vec2& v) const override;
	Vec3 Sample_Li(const Intersection& isect, Vec3* wi, real* pdf, Intersection* lightPoint, const Vec2& u) const override;
	real Pdf_Li(const Intersection& isect, const Vec3& wi) const override;
	Vec3 Power() const override;
	bool IsEnvironmentLight() const override { return true; }

private:
	std::unique_ptr<EnvironmentMap> envMap;
	Vec3 center;
	real radius;
};

NAMESPACE_END