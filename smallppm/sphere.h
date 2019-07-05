#pragma once

#include "utils.h"
#include "linagl.h"
#include "shape.h"
#include "ray.h"
#include "sampling.h"
#include "geometry_util.h"

NAMESPACE_BEGIN

class Sphere : public Shape {
public:
	Sphere(real radius, Vec3 position, Vec3 emission, Vec3 color, ReflectionType reflType) :
		rad(radius), p(position), Shape(reflType, color, emission, emission != Vec3()) {}

	real Intersect(const Ray &r, Intersection *isect) const override {
		// ray-sphere Intersection returns distance
		Vec3 op = p - r.o;
		real t, b = op.Dot(r.d), det = b * b - op.Dot(op) + rad * rad;
		//float64 b = (float64)op.x * (float64)r.d.x + (float64)op.y * (float64)r.d.y + (float64)op.z * (float64)r.d.z;
		//float64 det = b * b - ((float64)op.x * (float64)op.x + (float64)op.y * (float64)op.y + (float64)op.z * (float64)op.z)
		//	+ (float64)rad * (float64)rad;
		//real t;
		if (det < 0) {
			return r.tMax;
		}
		else {
			det = sqrt(det);
		}
		t = (t = b - det) > 1e-4 ? t : ((t = b + det) > 1e-4 ? t : r.tMax);

		isect->hit = r.o + r.d * t;
		isect->n = (isect->hit - p).Norm();
		isect->nl = isect->n.Dot(r.d) < 0 ? isect->n : isect->n * -1;
		isect->wo = -1 * r.d;
		isect->rayEps = 5e-4f * t;
		return t;
	}

	bool Intersect(const Ray &r, Intersection *isect, real *t) const override {
		// ray-sphere Intersection returns distance
		Vec3 op = p - r.o;
		real b = op.Dot(r.d), det = b * b - op.Dot(op) + rad * rad;
		//float64 b = (float64)op.x * (float64)r.d.x + (float64)op.y * (float64)r.d.y + (float64)op.z * (float64)r.d.z;
		//float64 det = b * b - ((float64)op.x * (float64)op.x + (float64)op.y * (float64)op.y + (float64)op.z * (float64)op.z) 
		//	+ (float64)rad * (float64)rad;
		//float64 tt;
		if (det < 0) {
			*t = r.tMax;
			return false;
		}
		else {
			det = std::sqrt(det);
		}

		*t = ((*t) = b - det) > 1e-6 ? (*t) : (((*t) = b + det) > 1e-6 ? (*t) : r.tMax);
		//tt = ((tt) = b - det) > 1e-6 ? (tt) : (((tt) = b + det) > 1e-6 ? (tt) : r.tMax);
		//*t = tt;
		//isect->hit = r.o + r.d * tt;
		isect->hit = r.o + r.d * (*t);
		isect->n = (isect->hit - p).Norm();
		isect->nl = isect->n.Dot(r.d) < 0 ? isect->n : isect->n * -1;
		isect->wo = -1 * r.d;
		isect->rayEps = 5e-4f * (*t);

		return (*t > r.tMin && *t < r.tMax);
	}


	bool Intersect(const Ray &r) const override {
		// ray-sphere Intersection returns distance
		Vec3 op = p - r.o;
		real t, b = op.Dot(r.d), det = b * b - op.Dot(op) + rad * rad;
		//float64 b = (float64)op.x * (float64)r.d.x + (float64)op.y * (float64)r.d.y + (float64)op.z * (float64)r.d.z;
		//float64 det = b * b - ((float64)op.x * (float64)op.x + (float64)op.y * (float64)op.y + (float64)op.z * (float64)op.z)
		//	+ (float64)rad * (float64)rad;
		//float64 t;
		if (det < 0) {
			return false;
		}
		else {
			det = sqrt(det);
		}
		t = (t = b - det) > 1e-4 ? t : ((t = b + det) > 1e-4 ? t : r.tMax);

		return t > r.tMin  && t < r.tMax;
	}

	Vec3 Sample(real *pdf, const Vec2 &u) const override {
		*pdf = 1.f / (4.f * PI * rad * rad);
		return UniformSampleSphere(u) * rad + p;
	}

	Vec3 Sample(const Intersection &isect, real *pdf, const Vec2 &u) const override {
		/*
		Vec3 sw = p - isect.hit, su = ((fabs(sw.x) > .1 ? Vec3(0, 1) : Vec3(1)) % sw).Norm(), sv = sw % su;
		real cos_a_max = sqrt(1 - rad * rad / (isect.hit - p).Dot(isect.hit - p));
		real zeta1 = u.x, zeta2 = u.y;
		real cos_a = 1 - zeta1 + zeta1 * cos_a_max;
		real sin_a = sqrt(1 - cos_a * cos_a);
		real phi = 2 * PI * zeta2;
		Vec3 dir = su * cos(phi) * sin_a + sv * sin(phi) * sin_a + sw * cos_a;
		real omega = 2 * PI *(1 - cos_a_max);
		*pdf = 1.0 / omega;
		return isect.hit + dir.Norm() * (sw.Length() - rad);
		//return dir.Norm();*/
		if ((p - isect.hit).Length2() <= rad * rad) {
			Vec3 lightPoint = Sample(pdf, u);
			Vec3 wi = lightPoint - isect.hit;
			if (wi.Dot(wi) == 0)
				*pdf = 0;
			else {
				real s = wi.Length();
				wi.Normalize();
				*pdf *= s / std::abs((lightPoint - p).Norm().Dot(-1 * wi));
			}
			if (std::isinf(*pdf)) *pdf = 0.f;
			return lightPoint;
		}
		Vec3 localZ = (p - isect.hit);
		real dis = localZ.Length();
		localZ.Normalize();
		Vec3 localX, localY;
		CoordinateSystem(localZ, &localX, &localY);
		real sin2ThetaMax = rad * rad / (dis * dis);
		real cosThetaMax = std::sqrt(std::max(1 - sin2ThetaMax, (real)0));
		real cosTheta = (1 - u[0]) + u[0] * cosThetaMax;
		real sinTheta = std::sqrt(std::max(1 - cosTheta * cosTheta, (real)0));
		real phi = 2 * PI * u[1];
		real s = dis * cosTheta - std::sqrt(std::max(rad * rad - dis * dis * sinTheta * sinTheta, (real)0));
		real cosAlpha = (dis * dis + rad * rad - s * s) / (2 * dis * rad);
		real sinAlpha = std::sqrt(std::max(1 - cosAlpha * cosAlpha, (real)0));
		Vec3 wi = sinAlpha * std::cos(phi) * localX + sinAlpha * std::sin(phi) * localY + cosTheta * localZ;
		Vec3 nWorld = -1 * sinAlpha * std::cos(phi) * localX - sinAlpha * std::sin(phi) * localY - cosAlpha * localZ;
		Vec3 lightPoint = p + rad * nWorld + nWorld * rayeps;
		//Vec3 wi = sinTheta * std::cos(phi) * localX + sinTheta * std::sin(phi) * localY + cosTheta * localZ;
		//Vec3 lightPoint = isect.hit + wi * s;
		*pdf = 1 / (2 * PI * (1 - cosThetaMax));

		return lightPoint;
	}

	Vec3 GetNorm(const Vec3 & point) const override {
		return (point - p).Norm();
	}

	real Area() const override {
		return 4.f * PI * rad * rad;
	}

private:
	real rad; Vec3 p;
};

NAMESPACE_END