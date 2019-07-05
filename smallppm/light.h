#pragma once

#include "linagl.h"
#include "intersection.h"
#include "shape.h"

NAMESPACE_BEGIN

class Light {
public:
	Light() {}
	virtual Vec3 DirectIllumination(const Intersection &isect, const std::shared_ptr<BSDF> &bsdf,
		const Vec3 &importance, Vec3 *dir, Intersection *lightPoint, const Vec2 &u) const = 0;
	virtual Vec3 Emission() const = 0;
	//virtual Vec3 SampleLight(Vec3 *pos, Vec3 *dir, Vec3 *lightNorm, real *pdfPos, real *pdfDir, const Vec2 &u, const Vec2 &v) const {
	//	SampleOnLight(pos, dir, lightNorm, pdfPos, pdfDir, u, v);
	//	return Emission();
	//}
	virtual Vec3 SampleLight(Intersection *isect, Vec3 *dir, real *pdfPos, real *pdfDir, const Vec2 &u, const Vec2 &v) const {
		SampleOnLight(isect, dir, pdfPos, pdfDir, u, v);
		return Emission();
	}
	virtual Vec3 Sample_Li(const Intersection &isect, Vec3 *wi, real *pdf, Intersection *lightPoint, const Vec2 &u) const = 0;
	virtual real Pdf_Li(const Intersection &isect, const Vec3 &wi) const = 0;
	virtual int GetId() const = 0;
	virtual Vec3 Power() const = 0;
	virtual bool IsAreaLight() const { return false; }
	virtual bool IsDeltaLight() const { return false; }
	virtual std::shared_ptr<Shape> GetShapePtr() const = 0;
protected:
	//virtual void SampleOnLight(Vec3 *pos, Vec3 *dir, Vec3 *lightNorm, real *pdfPos, real *pdfDir, const Vec2 &u, const Vec2 &v) const = 0;
	virtual void SampleOnLight(Intersection *isect, Vec3 *dir, real *pdfPos, real *pdfDir, const Vec2 &u, const Vec2 &v) const = 0;

};

NAMESPACE_END