#pragma once
#include "def.h"
#include "shape.h"
#include "transform.h"

NAMESPACE_BEGIN

class Triangle : public Shape {
public:
	Triangle(){}
	Triangle(Transform* ObjectToWorld, Transform* WorldToObject,
		const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& normal, Vec2 *uvs = nullptr) : 
		Shape(ObjectToWorld, WorldToObject), p0(p0), p1(p1), p2(p2), faceNormal(normal)
	{
		//faceNormal = Cross(p1 - p0, p2 - p0).Norm();
		n0 = n1 = n2 = faceNormal;
		if (uvs) {
			this->uvs[0] = uvs[0];
			this->uvs[1] = uvs[1];
			this->uvs[2] = uvs[2];
		}
		Initialize();

	}
	Triangle(Transform* ObjectToWorld, Transform* WorldToObject,
		const Vec3& p0, const Vec3& p1, const Vec3& p2, 
		const Vec3& n0, const Vec3& n1, const Vec3& n2, 
		Vec2* uvs = nullptr) :
		Shape(ObjectToWorld, WorldToObject), 
		p0(p0), p1(p1), p2(p2), n0(n0), n1(n1), n2(n2)
	{
		faceNormal = Cross(p1 - p0, p2 - p0).Norm();
		if (uvs) {
			this->uvs[0] = uvs[0];
			this->uvs[1] = uvs[1];
			this->uvs[2] = uvs[2];
		}
		Initialize();
	}

	Triangle(Transform* ObjectToWorld, Transform* WorldToObject,
		const Vec3& p0, const Vec3& p1, const Vec3& p2,
		const Vec3& n0, const Vec3& n1, const Vec3& n2,
		const Vec2 &uv0, const Vec2& uv1, const Vec2& uv2) :
		Shape(ObjectToWorld, WorldToObject),
		p0(p0), p1(p1), p2(p2), n0(n0), n1(n1), n2(n2)
	{
		faceNormal = Cross(p1 - p0, p2 - p0).Norm();
		this->uvs[0] = uv0;
		this->uvs[1] = uv1;
		this->uvs[2] = uv2;
		Initialize();
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

	AABB ObjectBound() const override {
		AABB ret(p0, p1);
		return Union(ret, p2);
	}

	void SetTransform(const Transform& transform) {
		p0 = transform(p0);
		p1 = transform(p1);
		p2 = transform(p2);
		n0 = transform.TransformNormal(n0);
		n1 = transform.TransformNormal(n1);
		n2 = transform.TransformNormal(n2);
		faceNormal = transform.TransformNormal(faceNormal);

		Initialize();
	}

	void QueryIntersectionInfo(const Ray& ray, Intersection* isct) const override;

private:
	void GetUVs(real uv[3][2]) const {
		if (uvs == nullptr) {
			uv[0][0] = 0.; uv[0][1] = 0.;
			uv[1][0] = 1.; uv[1][1] = 0.;
			uv[2][0] = 1.; uv[2][1] = 1.;
		}
		else {
			uv[0][0] = uvs[0][0]; uv[0][1] = uvs[0][1];
			uv[1][0] = uvs[1][0]; uv[1][1] = uvs[1][1];
			uv[2][0] = uvs[2][0]; uv[2][1] = uvs[2][1];
		}

	}

	void Initialize() {
		e1 = p1 - p0;
		e2 = p2 - p0;
		real uvs[3][2];
		GetUVs(uvs);
		du1 = uvs[0][0] - uvs[2][0];
		du2 = uvs[1][0] - uvs[2][0];
		dv1 = uvs[0][1] - uvs[2][1];
		dv2 = uvs[1][1] - uvs[2][1];
		dp1 = p0 - p2;
		dp2 = p1 - p2;
	}

	Vec3 p0, p1, p2;
	Vec3 n0, n1, n2;
	Vec3 faceNormal;
	Vec2 uvs[3];
	Vec3 e1, e2;

private:
	real du1; 
	real du2;
	real dv1;
	real dv2;
	Vec3 dp1;
	Vec3 dp2;
};

NAMESPACE_END
