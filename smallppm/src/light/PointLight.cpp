#include "PointLight.h"
#include "visual/Sampling.h"

NAMESPACE_BEGIN

Vec3 PointLight::Sample_Li(const Intersection& isect, Vec3* wi, real* pdf, Intersection* lightPoint, const Vec2& u) const {
	lightPoint->hit = pLight;
	Vec3 dir = isect.hit - pLight;
	*wi = -dir.Norm();
	lightPoint->n = lightPoint->nl = lightPoint->ng = -(*wi);
	*pdf = 1.f;
	return I / dir.Length2();
}


Vec3 PointLight::SampleLight(Intersection* isect, Vec3* dir, real* pdfPos, real* pdfDir, const Vec2& u, const Vec2& v) const {
	*pdfPos = 1;
	Vec3 sampledDir = UniformSampleSphere(v);
	*dir = sampledDir.Norm();
	*pdfDir = UniformSpherePdf();
	isect->hit = pLight;
	isect->n = isect->ng = isect->nl = *dir;
	return I;
}

NAMESPACE_END