#include "EnvLight.h"
#include "visual/Scene.h"

GY_NAMESPACE_BEGIN

void EnvironmentLight::Initialize(const Scene& scene) 
{
	scene.GetBoundingSphere(&mCenter, &mRadius);
}

Vec3 EnvironmentLight::Emission(const Ray& ray) const
{
	Vec3 w = ray.mDir.Norm();
	return mpEnvMap->Lookup(w);
}

Vec3 EnvironmentLight::Emission(const Vec3& dir, const Vec3 &normal) const 
{
	return mpEnvMap->Lookup(dir.Norm());
}

Vec3 EnvironmentLight::SampleLight(Intersection* isect, Vec3* dir, real* pdfPos,
	real* pdfDir, const Vec2& u, const Vec2& v) const 
{
	Vec3 sampledDir;
	Vec3 radiance = mpEnvMap->Sample(&sampledDir, pdfDir, u);
	*dir = -sampledDir;
	isect->n = *dir;

	Vec3 v1, v2;
	CoordinateSystem(-*dir, &v1, &v2);
	Vec2 cd = ConcentricSampleDisk(v);
	Vec3 pDisk = mCenter + mRadius * (cd.x * v1 + cd.y * v2);
	isect->hit = pDisk + mRadius * (-*dir);

	*pdfPos = 1 / (PI * mRadius * mRadius);
	//std::cout << radiance << " pdfDir: " << *pdfDir << " pdfPos: " << *pdfPos << std::endl;
	return radiance;
}


Vec3 EnvironmentLight::Sample_Li(const Intersection& isect, Vec3* wi, real* pdf,
	Intersection* lightPoint, const Vec2& u) const 
{
	Vec3 sampledDir;
	Vec3 radiance = mpEnvMap->Sample(&sampledDir, pdf, u);
	*wi = sampledDir.Norm();
	lightPoint->hit = isect.hit + *wi * (2 * mRadius);
	lightPoint->n = lightPoint->nl = lightPoint->ng = *wi;
	return radiance;
}

real EnvironmentLight::Pdf_Li(const Intersection& isect, const Vec3& wi) const 
{
	return mpEnvMap->PdfDir(wi);
}


Vec3 EnvironmentLight::Power() const 
{
	return mpEnvMap->LookupRadiance(0.5, 0.5) * PI * mRadius * mRadius;
}


GY_NAMESPACE_END