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
	Sphere(real radius, Vec position, Vec emission, Vec color, ReflectionType reflType) :
		rad(radius), p(position), Shape(reflType, color, emission, emission != Vec()) {}

	real Intersect(const Ray &r, Intersection *isect) const {
		// ray-sphere Intersection returns distance
		Vec op = p - r.o;
		real t, b = op.dot(r.d), det = b * b - op.dot(op) + rad * rad;
		if (det < 0) {
			return Inf;
		}
		else {
			det = sqrt(det);
		}
		t = (t = b - det) > 1e-4 ? t : ((t = b + det) > 1e-4 ? t : Inf);

		isect->hit = r.o + r.d * t;
		isect->n = (isect->hit - p).norm();
		isect->nl = isect->n.dot(r.d) < 0 ? isect->n : isect->n * -1;
		isect->wo = -1 * r.d;
		return t;
	}

	bool Intersect(const Ray &r, Intersection *isect, real *t) const {
		// ray-sphere Intersection returns distance
		Vec op = p - r.o;
		real b = op.dot(r.d), det = b * b - op.dot(op) + rad * rad;
		if (det < 0) {
			*t = Inf;
			return false;
		}
		else {
			det = sqrt(det);
		}
		*t = ((*t) = b - det) > 1e-4 ? (*t) : (((*t) = b + det) > 1e-4 ? (*t) : Inf);

		isect->hit = r.o + r.d * (*t);
		isect->n = (isect->hit - p).norm();
		isect->nl = isect->n.dot(r.d) < 0 ? isect->n : isect->n * -1;
		isect->wo = -1 * r.d;

		return (*t > 0 && *t < r.tMax);
	}


	bool Intersect(const Ray &r) const {
		// ray-sphere Intersection returns distance
		Vec op = p - r.o;
		real t, b = op.dot(r.d), det = b * b - op.dot(op) + rad * rad;
		if (det < 0) {
			return false;
		}
		else {
			det = sqrt(det);
		}
		t = (t = b - det) > 1e-4 ? t : ((t = b + det) > 1e-4 ? t : Inf);

		return t > 0 && t < r.tMax;
	}

	Vec Sample(real *pdf, Vec u) const {
		*pdf = 1.0 / (4.0 * PI * rad * rad);
		return UniformSampleSphere(u) * rad + p;
	}

	Vec Sample(const Intersection &isect, real *pdf, Vec u) const {
		/*
		Vec sw = p - isect.hit, su = ((fabs(sw.x) > .1 ? Vec(0, 1) : Vec(1)) % sw).norm(), sv = sw % su;
		real cos_a_max = sqrt(1 - rad * rad / (isect.hit - p).dot(isect.hit - p));
		real zeta1 = u.x, zeta2 = u.y;
		real cos_a = 1 - zeta1 + zeta1 * cos_a_max;
		real sin_a = sqrt(1 - cos_a * cos_a);
		real phi = 2 * PI * zeta2;
		Vec dir = su * cos(phi) * sin_a + sv * sin(phi) * sin_a + sw * cos_a;
		real omega = 2 * PI *(1 - cos_a_max);
		*pdf = 1.0 / omega;
		return isect.hit + dir.norm() * (sw.length() - rad);
		//return dir.norm();*/
		if ((p - isect.hit).length2() <= rad * rad) {
			Vec lightPoint = Sample(pdf, u);
			Vec wi = lightPoint - isect.hit;
			if (wi.dot(wi) == 0)
				*pdf = 0;
			else {
				real s = wi.length();
				wi.normalize();
				*pdf *= s / std::abs((lightPoint - p).norm().dot(-1 * wi));
			}
			if (std::isinf(*pdf)) *pdf = 0.f;
			return lightPoint;
		}
		Vec localZ = (p - isect.hit);
		real dis = localZ.length();
		localZ.normalize();
		Vec localX, localY;
		CoordinateSystem(localZ, &localX, &localY);
		real sin2ThetaMax = rad * rad / (dis * dis);
		real cosThetaMax = std::sqrt(std::max(1 - sin2ThetaMax, (real)0));
		real cosTheta = (1 - u[0]) + u[0] * cosThetaMax;
		real sinTheta = std::sqrt(std::max(1 - cosTheta * cosTheta, (real)0));
		real phi = 2 * PI * u[1];
		real s = dis * cosTheta - std::sqrt(std::max(rad * rad - dis * dis * sinTheta * sinTheta, (real)0));
		real cosAlpha = (dis * dis + rad * rad - s * s) / (2 * dis * rad);
		real sinAlpha = std::sqrt(std::max(1 - cosAlpha * cosAlpha, (real)0));
		Vec wi = sinAlpha * std::cos(phi) * localX + sinAlpha * std::sin(phi) * localY + cosTheta * localZ;
		Vec nWorld = -1 * sinAlpha * std::cos(phi) * localX - sinAlpha * std::sin(phi) * localY - cosAlpha * localZ;
		Vec lightPoint = p + rad * nWorld + nWorld * rayeps;
		//Vec wi = sinTheta * std::cos(phi) * localX + sinTheta * std::sin(phi) * localY + cosTheta * localZ;
		//Vec lightPoint = isect.hit + wi * s;
		*pdf = 1 / (2 * PI * (1 - cosThetaMax));

		return lightPoint;
	}

	Vec GetNorm(const Vec & point) const {
		return (point - p).norm();
	}

	real Area() const {
		return 4.0 * PI * rad * rad;
	}

private:
	real rad; Vec p;
};

NAMESPACE_END