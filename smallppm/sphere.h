#pragma once

#include "utils.h"
#include "linagl.h"
#include "shape.h"
#include "ray.h"
#include "sampling.h"
#include "geometry_util.h"
#include "debug_utils.h"

NAMESPACE_BEGIN

class Sphere : public Shape {
public:
	Sphere(real radius, Vec3 position) : Shape(), rad(radius), p(position) {}

	bool Intersect(const Ray &r, Intersection *isect, real *t) const override {
		// ray-sphere Intersection returns distance
		Vec3 v = r.o - p;
		real vd = v.Dot(r.d), d2 = r.d.Dot(r.d), v2 = v.Dot(v);
		real det = vd * vd - d2 * v2 + d2 * rad * rad;


		if (det < 0) {
			*t = r.tMax;
			return false;
		}
		else {
			det = std::sqrt(det);
		}

		*t = (-vd - det) / d2;
		if (*t < 1e-4) {
			*t = (-vd + det) / d2;
			if (*t < 1e-4) {
				*t = r.tMax;
			}
		}

		isect->hit = r.o + r.d * (*t);
		Vec3 scaledDir = (isect->hit - p) * rad / Distance(isect->hit, p);
		isect->hit = p + scaledDir;
		isect->n = (isect->hit - p).Norm();
		isect->nl = isect->n.Dot(r.d) < 0 ? isect->n : isect->n * -1;
		isect->wo = -1 * r.d;
		//isect->pError = gamma(5) * Abs(isect->hit);
		//isect->pError = Vec3(1e-4f, 1e-4f, 1e-4f);

		//real vdError = vd * gamma(2 + 1);
		//real sqrtDetError = vd * gamma(4) + d2 * v2 * gamma(4) + d2 * rad * rad * gamma(1);
		//real tError = vd / d2 * gamma(3) + sqrtDetError / d2 * gamma(2);
		//Vec3 hitError = gamma(1) * Abs(r.o) + tError * (1 + gamma(2)) * Abs(r.d) + gamma(2) * Abs(tt * r.d);
		//isect->pError = hitError;
		isect->pError = Abs(p) * gamma(1) + Abs(scaledDir) * gamma(6);
		//std::cout << isect->pError << std::endl;
		return (*t > r.tMin && *t < r.tMax);
	}


	bool Intersect(const Ray &r) const override {
		// ray-sphere Intersection returns distance
		Vec3 op = p - r.o;

		Vec3 v = r.o - p;
		real vd = v.Dot(r.d), d2 = r.d.Dot(r.d), v2 = v.Dot(v);
		real t;
		real det = vd * vd - d2 * v.Dot(v) + d2 * rad * rad;

		if (det < 0) {
			return false;
		}
		else {
			det = sqrt(det);
		}
		t = (-vd - det) / d2;
		if (t < 1e-4) {
			t = (-vd + det) / d2;
			if (t < 1e-4) {
				t = r.tMax;
			}
		}

		return (t > r.tMin  && t < r.tMax);
	}

	Intersection Sample(real *pdf, const Vec2 &u) const override {
		*pdf = 1.f / (4.f * PI * rad * rad);
		Vec3 pHit = UniformSampleSphere(u) * rad + p;
		Vec3 scaledDir((pHit - p) * rad / Distance(pHit, p));
		Intersection isect;
		isect.hit = p + scaledDir;
		isect.n = (pHit - p).Norm();
		isect.nl = isect.n;

		Vec3 pError = Abs(p) * gamma(1) + Abs(scaledDir) * gamma(6);
		isect.pError = pError;
		return isect;
	}

	Intersection Sample(const Intersection &isect, real *pdf, const Vec2 &u) const override {
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
			Intersection lightPoint = Sample(pdf, u);
			Vec3 wi = lightPoint.hit - isect.hit;
			if (wi.Dot(wi) == 0)
				*pdf = 0;
			else {
				real s = wi.Length();
				wi.Normalize();
				*pdf *= s / std::abs((lightPoint.hit - p).Norm().Dot(-1 * wi));
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
		Vec3 pWorld = p + rad * nWorld + nWorld * rayeps;
		//Vec3 wi = sinTheta * std::cos(phi) * localX + sinTheta * std::sin(phi) * localY + cosTheta * localZ;
		//Vec3 lightPoint = isect.hit + wi * s;
		*pdf = 1 / (2 * PI * (1 - cosThetaMax));
		Intersection ret;
		ret.hit = pWorld;
		ret.n = nWorld;
		ret.nl = ret.n;
		ret.pError = gamma(5) * Abs(pWorld);
		//std::cout << ret.hit << " " << ret.pError << std::endl;
		return ret;
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