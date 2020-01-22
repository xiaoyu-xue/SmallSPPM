#include "environment_light.h"
#include "scene.h"

NAMESPACE_BEGIN

void EnvironmentLight::Initialize(const Scene& scene) {
	scene.GetBoundingSphere(&center, &radius);
}

Vec3 EnvironmentLight::Emission(const Ray& ray) const {
	Vec3 w = ray.d.Norm();
	return envMap->Lookup(w);
}

Vec3 EnvironmentLight::Emission(const Vec3& dir, const Vec3 &normal) const {
	return envMap->Lookup(dir.Norm());
}

Vec3 EnvironmentLight::SampleLight(Intersection* isect, Vec3* dir, real* pdfPos,
	real* pdfDir, const Vec2& u, const Vec2& v) const {
	Vec3 sampledDir;
	Vec3 radiance = envMap->Sample(&sampledDir, pdfDir, u);
	*dir = -sampledDir;
	isect->n = *dir;

	Vec3 v1, v2;
	CoordinateSystem(-*dir, &v1, &v2);
	Vec2 cd = ConcentricSampleDisk(v);
	Vec3 pDisk = center + radius * (cd.x * v1 + cd.y * v2);
	isect->hit = pDisk + radius * (-*dir);

	*pdfPos = 1 / (PI * radius * radius);
	//std::cout << radiance << " pdfDir: " << *pdfDir << " pdfPos: " << *pdfPos << std::endl;
	return radiance;
}


Vec3 EnvironmentLight::Sample_Li(const Intersection& isect, Vec3* wi, real* pdf,
	Intersection* lightPoint, const Vec2& u) const {
	Vec3 sampledDir;
	Vec3 radiance = envMap->Sample(&sampledDir, pdf, u);
	*wi = sampledDir.Norm();
	lightPoint->hit = isect.hit + *wi * (2 * radius);
	lightPoint->n = lightPoint->nl = lightPoint->ng = *wi;
	return radiance;
}

real EnvironmentLight::Pdf_Li(const Intersection& isect, const Vec3& wi) const {
	return envMap->PdfDir(wi);
}


Vec3 EnvironmentLight::Power() const {
	return envMap->LookupRadiance(0.5, 0.5) * PI * radius * radius;
}


NAMESPACE_END