#pragma once

#include "linagl.h"
#include "intersection.h"
#include "shape.h"

class Light {
public:
	Light() {}
	virtual Vec DirectIllumination(const Intersection &isect, const std::shared_ptr<BSDF> &bsdf,
		const Vec &importance, Vec *dir, Intersection *lightPoint, Vec u) const = 0;
	virtual Vec Emission() const = 0;
	virtual Vec SampleLight(Vec *pos, Vec *dir, Vec *lightNorm, real *pdfPos, real *pdfDir, Vec u, Vec v) const {
		SampleOnLight(pos, dir, lightNorm, pdfPos, pdfDir, u, v);
		return Emission();
	}
	virtual Vec Sample_Li(const Intersection &isect, Vec *wi, real *pdf, Intersection *lightPoint, Vec u) const = 0;
	virtual real Pdf_Li(const Intersection &isect, const Vec &wi) const = 0;
	virtual int GetId() const = 0;
	virtual Vec Power() const = 0;
	virtual bool IsAreaLight() const { return false; }
	virtual bool IsDeltaLight() const { return false; }
	virtual std::shared_ptr<Shape> GetShapePtr() const = 0;
protected:
	virtual void SampleOnLight(Vec *pos, Vec *dir, Vec *lightNorm, real *pdfPos, real *pdfDir, Vec u, Vec v) const = 0;
};
