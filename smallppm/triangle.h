#pragma once
#include "def.h"
#include "shape.h"
#include "transform.h"

NAMESPACE_BEGIN

class Triangle : public Shape {
public:
	Triangle(Transform* ObjectToWorld, Transform* WorldToObject,
		const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& normal) : 
		Shape(ObjectToWorld, WorldToObject), p0(p0), p1(p1), p2(p2), faceNormal(normal)
	{

	}

	bool Intersect(const Ray& r, Intersection* isect, real* t) const override;

	bool Intersect(const Ray& r) const override;

	Intersection Sample(real* pdf, const Vec2& rand) const override;

	//Intersection Sample(const Intersection& isect, real* pdf, const Vec2& u) const override;

	Vec3 GetNorm(const Vec3& point) const override {
		return faceNormal;
	}

	real Area() const override {
		return 0.5 * Cross(p1 - p0, p2 - p0).Length();
	}
private:
	void GetUVs(real uv[3][2]) const {
		uv[0][0] = 0.; uv[0][1] = 0.;
		uv[1][0] = 1.; uv[1][1] = 0.;
		uv[2][0] = 1.; uv[2][1] = 1.;
	}

	Vec3 p0, p1, p2;
	Vec3 faceNormal;
};

NAMESPACE_END
