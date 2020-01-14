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

	Sphere(const Transform* ObjectToWorld, const Transform* WorldToObject,
		real radius, Vec3 position, bool inside = false) :
		Shape(ObjectToWorld, WorldToObject), rad(radius), p(position), inside(inside) {

		aabb = AABB(p - Vec3(rad, rad, rad), p + Vec3(rad, rad, rad));
	}

	bool Intersect(const Ray &ray, Intersection *isect, real *t) const override {
		// ray-sphere Intersection returns distance

		Ray r;
		real tHit;
		if (WorldToObject) {
			r = (*WorldToObject)(ray);
		}
		else {
			r = ray;
		}

		Vec3 po = r.o - p;
		real vd = po.Dot(r.d), d2 = r.d.Dot(r.d), po2 = po.Dot(po);
		real det = vd * vd - d2 * po2 + d2 * rad * rad;


		if (det < 0) {
			//*t = r.tMax;
			return false;
		}
		else {
			det = std::sqrt(det);
		}

		tHit = (-vd - det) / d2;
		if (tHit < r.tMin || tHit > r.tMax) {
			tHit = (-vd + det) / d2;
			if (tHit < r.tMin || tHit > r.tMax) {
				return false;
			}
		}

		//DEBUG_PIXEL_IF() {
		//	std::cout << "Intersection function: tmax: " << r.tMax << " tHit: " << tHit << std::endl;
		//}

		//isect->hit = r.o + r.d * (*t);
		//Vec3 scaledDir = (isect->hit - p) * rad / Distance(isect->hit, p);
		//isect->hit = p + scaledDir;
		//isect->n = (isect->hit - p).Norm();
		//isect->nl = isect->n.Dot(r.d) < 0 ? isect->n : isect->n * -1;
		//isect->wo = -1 * r.d;
		//isect->pError = Abs(p) * gamma(1) + Abs(scaledDir) * gamma(6);

		Vec3 hit = r.o + r.d * (tHit);
		Vec3 scaledDir = (hit - p) * rad / Distance(hit, p);
		hit = p + scaledDir;
		Vec3 n = GetNorm(hit);
		Vec3 nl = n.Dot(r.d) < 0 ? n : n * -1;
		Vec3 wo = -1 * r.d;
		Vec3 pError = Abs(p) * gamma(1) + Abs(scaledDir) * gamma(6);

		//Compute dpdu, dpdv
		Vec3 pHit = hit - p;
		real phi = std::atan2(pHit.y, pHit.x);
		if (phi < 0) phi += 2 * PI;
		real u = phi / phiMax;
		real theta = std::acos(Clamp(pHit.z / rad, -1, 1));
		real v = (theta - thetaMin) / (thetaMax - thetaMin);

		real zRadius = std::sqrt(pHit.x * pHit.x + pHit.y * pHit.y);
		real invZRadius = 1 / zRadius;
		real cosPhi = pHit.x * invZRadius;
		real sinPhi = pHit.y * invZRadius;
		Vec3 dpdu(-phiMax * pHit.y, phiMax * pHit.x, 0);
		Vec3 dpdv = (thetaMax - thetaMin) * Vec3(pHit.z * cosPhi, pHit.z * sinPhi, -rad * std::sin(theta));
		//CoordinateSystem(GetNorm(hit), &dpdu, &dpdv);

		Intersection it(hit, n, nl, dpdu, dpdv, wo, pError);
		it.SetShading(n, dpdu, dpdv);
		if (ObjectToWorld) {
			*isect = (*ObjectToWorld)(it);
		}
		else {
			*isect = it;
		}
		isect->uv = Vec2(u, v);
		*t = tHit;

		return true;

	}


	bool Intersect(const Ray &ray) const override {
		// ray-sphere Intersection returns distance

		Ray r;
		real tHit;
		if (WorldToObject) {
			r = (*WorldToObject)(ray);
		}
		else {
			r = ray;
		}

		Vec3 op = p - r.o;
		Vec3 v = r.o - p;
		real vd = v.Dot(r.d), d2 = r.d.Dot(r.d), v2 = v.Dot(v);
		real det = vd * vd - d2 * v.Dot(v) + d2 * rad * rad;

		if (det < 0) {
			return false;
		}
		else {
			det = sqrt(det);
		}
		tHit = (-vd - det) / d2;
		if (tHit < r.tMin || tHit > r.tMax) {
			tHit = (-vd + det) / d2;
			if (tHit < r.tMin || tHit > r.tMax) {
				return false;
			}
		}

		return true;
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
		Vec3 n = GetNorm(hit);
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
		isect.shapeId = shapeId;
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
		ret.shapeId = shapeId;
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
		//return (point - p).Norm();
		if (inside) return (p - point).Norm();
		else return (point - p).Norm();
	}

	real Area() const override {
		return 4.f * PI * rad * rad;
	}

	AABB ObjectBound() const override {
		return aabb;
	}

private:
	real rad; 
	Vec3 p;
	bool inside;
	AABB aabb;
	const real thetaMin = 0;
	const real thetaMax = PI;
	const real phiMax = 2 * PI;
};

NAMESPACE_END