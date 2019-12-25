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

	Sphere(real radius, Vec3 position, 
		const Transform* ObjectToWorld = nullptr, const Transform* WorldToObject = nullptr) :
		Shape(ObjectToWorld, WorldToObject), rad(radius), p(position) {

	}

	bool Intersect(const Ray &ray, Intersection *isect, real *t) const override {
		// ray-sphere Intersection returns distance

		Ray r;
		if (WorldToObject) {
			r = (*WorldToObject)(ray);
		}
		else {
			r = ray;
		}

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

		//isect->hit = r.o + r.d * (*t);
		//Vec3 scaledDir = (isect->hit - p) * rad / Distance(isect->hit, p);
		//isect->hit = p + scaledDir;
		//isect->n = (isect->hit - p).Norm();
		//isect->nl = isect->n.Dot(r.d) < 0 ? isect->n : isect->n * -1;
		//isect->wo = -1 * r.d;
		//isect->pError = Abs(p) * gamma(1) + Abs(scaledDir) * gamma(6);


		Vec3 hit = r.o + r.d * (*t);
		Vec3 scaledDir = (hit - p) * rad / Distance(hit, p);
		hit = p + scaledDir;
		Vec3 n = (hit - p).Norm();
		Vec3 nl = n.Dot(r.d) < 0 ? n : n * -1;
		Vec3 wo = -1 * r.d;
		Vec3 pError = Abs(p) * gamma(1) + Abs(scaledDir) * gamma(6);

		if (ObjectToWorld) {
			if (shapeId == 1) {
				//std::cout << "r.tMin: " << r.tMin << ", " << "t: "<< *t << ", " << "t.tMax: "<< r.tMax << std::endl;
				//std::cout << (*t > r.tMin&&* t < r.tMax) << std::endl;
			}
			*isect = (*ObjectToWorld)(Intersection(hit, n, nl, wo, pError));
			//if(shapeId == 0) std::cout << pError << std::endl << isect->pError << std::endl << std::endl;
			const Transform o2w = *ObjectToWorld;
			//*isect = Intersection(o2w(hit), o2w.TransformNormal(n), o2w.TransformNormal(nl), o2w.TransformVector(wo), pError);
		}
		else {
			*isect = Intersection(hit, n, nl, wo, pError);
		}


		//std::cout << isect->pError << std::endl;
		return (*t > r.tMin && *t < r.tMax);
	}


	bool Intersect(const Ray &ray) const override {
		// ray-sphere Intersection returns distance

		Ray r;
		if (WorldToObject) {
			r = (*WorldToObject)(ray);
		}
		else {
			r = ray;
		}

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
		//*pdf = 1.f / (4.f * PI * rad * rad);
		//Vec3 pHit = UniformSampleSphere(u) * rad + p;
		//Vec3 scaledDir((pHit - p) * rad / Distance(pHit, p));
		//Intersection isect;
		//isect.hit = p + scaledDir;
		//isect.n = (pHit - p).Norm();
		//isect.nl = isect.n;

		//Vec3 pError = Abs(p) * gamma(1) + Abs(scaledDir) * gamma(6);
		//isect.pError = pError;
		//return isect;



		*pdf = 1.f / (4.f * PI * rad * rad);
		Vec3 pHit = UniformSampleSphere(u) * rad + p;
		Vec3 scaledDir((pHit - p) * rad / Distance(pHit, p));
		Intersection isect;
		Vec3 hit = p + scaledDir;
		Vec3 n = (pHit - p).Norm();
		Vec3 nl = isect.n;
		Vec3 pError = Abs(p) * gamma(1) + Abs(scaledDir) * gamma(6);
		if (ObjectToWorld) {
			isect.hit = (*ObjectToWorld)(hit, pError, &isect.pError);
			isect.n = (*ObjectToWorld).TransformNormal(n).Norm();
		}
		else {
			isect.hit = hit;
			isect.n = n;
		}
		isect.nl = nl;
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


		//if ((p - isect.hit).Length2() <= rad * rad) {
		//	Intersection lightPoint = Sample(pdf, u);
		//	Vec3 wi = lightPoint.hit - isect.hit;
		//	if (wi.Dot(wi) == 0)
		//		*pdf = 0;
		//	else {
		//		real s = wi.Length();
		//		wi.Normalize();
		//		*pdf *= s / std::abs((lightPoint.hit - p).Norm().Dot(-1 * wi));
		//	}
		//	if (std::isinf(*pdf)) *pdf = 0.f;
		//	return lightPoint;
		//}
		//Vec3 localZ = (p - isect.hit);
		//real dis = localZ.Length();
		//localZ.Normalize();
		//Vec3 localX, localY;
		//CoordinateSystem(localZ, &localX, &localY);
		//real sin2ThetaMax = rad * rad / (dis * dis);
		//real cosThetaMax = std::sqrt(std::max(1 - sin2ThetaMax, (real)0));
		//real cosTheta = (1 - u[0]) + u[0] * cosThetaMax;
		//real sinTheta = std::sqrt(std::max(1 - cosTheta * cosTheta, (real)0));
		//real phi = 2 * PI * u[1];
		//real s = dis * cosTheta - std::sqrt(std::max(rad * rad - dis * dis * sinTheta * sinTheta, (real)0));
		//real cosAlpha = (dis * dis + rad * rad - s * s) / (2 * dis * rad);
		//real sinAlpha = std::sqrt(std::max(1 - cosAlpha * cosAlpha, (real)0));
		//Vec3 wi = sinAlpha * std::cos(phi) * localX + sinAlpha * std::sin(phi) * localY + cosTheta * localZ;
		//Vec3 nWorld = -1 * sinAlpha * std::cos(phi) * localX - sinAlpha * std::sin(phi) * localY - cosAlpha * localZ;
		//Vec3 pWorld = p + rad * nWorld + nWorld * rayeps;
		////Projection into surface and more convenient to calculate error
		//Vec3 scaledDir = (pWorld - pCenter) * rad / Distance(pWorld, pCenter);
		//pWorld = pWorld + scaledDir;
		////Vec3 wi = sinTheta * std::cos(phi) * localX + sinTheta * std::sin(phi) * localY + cosTheta * localZ;
		////Vec3 lightPoint = isect.hit + wi * s;
		//*pdf = 1 / (2 * PI * (1 - cosThetaMax));
		//Intersection ret;
		//ret.hit = pWorld;
		//ret.n = nWorld;
		//ret.nl = ret.n;
		//ret.pError = gamma(5) * Abs(pWorld);
		
		Vec3 pCenter;
		if (ObjectToWorld) {
			pCenter = (*ObjectToWorld)(p);
		}
		else {
			pCenter = p;
		}
		if ((pCenter - isect.hit).Length2() <= rad * rad) {
			Intersection lightPoint = Sample(pdf, u);
			Vec3 wi = lightPoint.hit - isect.hit;
			if (wi.Dot(wi) == 0)
				*pdf = 0;
			else {
				real s = wi.Length();
				wi.Normalize();
				*pdf *= s / std::abs((lightPoint.hit - pCenter).Norm().Dot(-1 * wi));
			}
			if (std::isinf(*pdf)) *pdf = 0.f;
			return lightPoint;
		}
		Vec3 localZ = (pCenter - isect.hit);
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
		nWorld.Normalize();
		Vec3 pWorld = pCenter + rad * nWorld / nWorld.Length() + rad * 1.5e-3 * nWorld.Norm();
		//Projection into surface and more convenient to calculate error
		//Vec3 scaledDir = (pWorld - pCenter) * rad / Distance(pWorld, pCenter);
		//pWorld = pWorld + scaledDir;
		//Vec3 wi = sinTheta * std::cos(phi) * localX + sinTheta * std::sin(phi) * localY + cosTheta * localZ;
		//Vec3 lightPoint = isect.hit + wi * s;
		*pdf = 1 / (2 * PI * (1 - cosThetaMax));
		Intersection ret;
		ret.hit = pWorld;
		ret.n = nWorld;
		ret.nl = ret.n;
		//ret.pError = gamma(6) * Abs(pWorld) + gamma(6) * Abs(scaledDir);
		ret.pError = gamma(5) * Abs(pWorld);
		return ret;


		//if ((p - isect.hit).Length2() <= rad * rad) {
		//	Intersection lightPoint = Sample(pdf, u);
		//	Vec3 wi = lightPoint.hit - isect.hit;
		//	if (wi.Dot(wi) == 0)
		//		*pdf = 0;
		//	else {
		//		real s = wi.Length();
		//		wi.Normalize();
		//		*pdf *= s / std::abs((lightPoint.hit - p).Norm().Dot(-1 * wi));
		//	}
		//	if (std::isinf(*pdf)) *pdf = 0.f;
		//	return lightPoint;
		//}
		//Vec3 localZ = (p - isect.hit);
		//real dis = localZ.Length();
		//localZ.Normalize();
		//Vec3 localX, localY;
		//CoordinateSystem(localZ, &localX, &localY);
		//real sin2ThetaMax = rad * rad / (dis * dis);
		//real cosThetaMax = std::sqrt(std::max(1 - sin2ThetaMax, (real)0));
		//real cosTheta = (1 - u[0]) + u[0] * cosThetaMax;
		//real sinTheta = std::sqrt(std::max(1 - cosTheta * cosTheta, (real)0));
		//real phi = 2 * PI * u[1];
		//real s = dis * cosTheta - std::sqrt(std::max(rad * rad - dis * dis * sinTheta * sinTheta, (real)0));
		//real cosAlpha = (dis * dis + rad * rad - s * s) / (2 * dis * rad);
		//real sinAlpha = std::sqrt(std::max(1 - cosAlpha * cosAlpha, (real)0));
		//Vec3 wi = sinAlpha * std::cos(phi) * localX + sinAlpha * std::sin(phi) * localY + cosTheta * localZ;
		//Vec3 nWorld = -1 * sinAlpha * std::cos(phi) * localX - sinAlpha * std::sin(phi) * localY - cosAlpha * localZ;
		//Vec3 pWorld = p + rad * nWorld + nWorld * rayeps;
		////Vec3 wi = sinTheta * std::cos(phi) * localX + sinTheta * std::sin(phi) * localY + cosTheta * localZ;
		////Vec3 lightPoint = isect.hit + wi * s;
		//*pdf = 1 / (2 * PI * (1 - cosThetaMax));
		//Intersection ret;
		//ret.hit = pWorld;
		//ret.n = nWorld;
		//ret.nl = ret.n;
		//ret.pError = gamma(5) * Abs(pWorld);
		////std::cout << ret.hit << " " << ret.pError << std::endl;
		//return ret;
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