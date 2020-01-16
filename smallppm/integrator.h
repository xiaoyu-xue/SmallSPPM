#pragma once

#include "utils.h"
#include "scene.h"
#include "intersection.h"
#include "visibility.h"

#include "debug_utils.h"

NAMESPACE_BEGIN

FORCE_INLINE real PowerHeuristic(int nf, real fPdf, int ng, real gPdf) {
	real f = nf * fPdf, g = ng * gPdf;
	return (f * f) / (f * f + g * g);
}

class Integrator {
public:
	virtual void Render(const Scene &scene) = 0;
	virtual ~Integrator() {}
	
	//static Vec3 DirectIllumination(const Scene &scene, const Intersection &isect, const std::shared_ptr<BSDF> &bsdf,
	//	real uLight, const Vec2 &u, const Vec3 &v) {
	//	Vec3 L;
	//	/*
	//	const std::vector<std::shared_ptr<Light>> lights = scene.GetLights();
	//	for (auto light : lights) {
	//		Vec3 dir;
	//		std::shared_ptr<Shape> hitObj;
	//		real t;
	//		Intersection lightPoint;
	//		Vec3 Li = light->DirectIllumination(isect, bsdf, importance, &dir, &lightPoint, v);
	//		Intersection intersection;
	//		if (scene.Intersect(Ray(isect.hit, dir), &t, &intersection, hitObj) && hitObj->GetId() == light->GetId()) {
	//			L = L + Li;
	//		}
	//	}
	//	return L;
	//	*/

	//	//Sample light
	//	real lightSamplingPdf;
	//	real weight1 = 1, weight2 = 1;
	//	Vec3 L1, L2;
	//	std::shared_ptr<Light> light = scene.SampleOneLight(&lightSamplingPdf, uLight);
	//	{
	//		Vec3 wi;
	//		Intersection lightPoint;
	//		real pdf;
	//		Vec3 Li = light->Sample_Li(isect, &wi, &pdf, &lightPoint, u);
	//		//{
	//		//	if (debugPixel == 1) {
	//		//		std::cout << "lightSamplePdf: " << pdf << std::endl;
	//		//	}
	//		//}
	//		if (pdf > 0) {
	//			real scatteringPdf = bsdf->Pdf(isect.wo, wi);
	//			Vec3 f = bsdf->f(isect.wo, wi);
	//			VisibilityTester visibilityTester(isect, lightPoint);
	//			//{ 
	//			//	{
	//			//		if (debugPixel == 1) {
	//			//			std::cout << Distance(lightPoint.hit, Vec3(50, 81.6f - 16.5f, 81.6f)) << std::endl;
	//			//			std::cout << "test n: " << isect.n.Dot(lightPoint.hit - isect.hit) 
	//			//				<< " " << isect.nl.Dot(lightPoint.hit - isect.hit) << std::endl;
	//			//		}
	//			//	}

	//			//}
	//			if (!visibilityTester.Unoccluded(scene)) {
	//				//if (debugPixel == 1) {
	//				//	//std::cout << "ok" << std::endl;
	//				//}
	//				weight1 = PowerHeuristic(1, pdf, 1, scatteringPdf);
	//				L1 = Li * f * std::abs(isect.n.Dot(wi)) / pdf;

	//				//if (debugPixel == 1) {
	//				//	std::cout << "sample light: weight " << weight1 << " cos: " << std::abs(isect.n.Dot(wi)) << " Li " << Li
	//				//		<< " pdf: " << pdf << std::endl;
	//				//}

	//				//std::cout << "ok" << std::endl;
	//				//std::cout << Li << std::endl;
	//				//debugPixel++;
	//			}
	//			//else {
	//			//	if (debugPixel == 1) {
	//			//		std::cout << "hitpoint: " << isect.hit << " " << isect.wo << std::endl;
	//			//		std::cout << "lightpoint: " << lightPoint.hit << std::endl << std::endl;
	//			//	}
	//			//}
	//		}

	//	}


	//	//Sample BSDF
	//	{
	//		Vec3 wi;
	//		real pdf;
	//		Vec3 f = bsdf->Sample_f(isect.wo, &wi, &pdf, v);
	//		Intersection intersection;
	//		std::shared_ptr<Shape> hitObj;
	//		real t = Inf;
	//		if (pdf != 0) {
	//			real lightPdf = light->Pdf_Li(isect, wi);

	//			//{
	//			//	if (debugPixel == 1) {
	//			//		std::cout << "lightSamplingPdf2: " << lightPdf << std::endl;

	//			//	}
	//			//}
	//			weight2 = PowerHeuristic(1, pdf, 1, lightPdf);
	//			if (scene.Intersect(isect.SpawnRay(wi), &intersection, &t) && isect.primitive->IsLight()) {
	//				std::shared_ptr<Light> emissionShape = isect.primitive->GetLight();
	//				L2 = emissionShape->Emission(intersection, -wi) * f * std::abs(isect.n.Dot(wi)) / pdf;
	//			}
	//			//if (scene.Intersect(Ray(isect.hit, wi, Inf, isect.rayEps), &t, &intersection, hitObj) && hitObj->IsLight()) {
	//			//	L2 = hitObj->GetEmission() * f * std::abs(isect.n.Dot(wi)) / pdf;
	//			//}
	//			/*if (scene.Intersect(Ray(isect.hit + wi * rayeps, wi), &t, &intersection, hitObj) && hitObj->IsLight()) {
	//				L2 = hitObj->GetEmission() * f * std::abs(isect.n.Dot(wi)) / pdf;
	//			}*/
	//		}
	//	}
	//	//return L1 / lightSamplingPdf;
	//	//return L2 / lightSamplingPdf;
	//	return (L1 * weight1 + L2 * weight2) / lightSamplingPdf;

	//}


	static Vec3 DirectIllumination(const Scene& scene, const Intersection& isect, BSDF* bsdf,
		real uLight, const Vec2& u, const Vec3& v) {
		Vec3 L;
		//Sample light
		real lightSamplingPdf;
		real weight1 = 1, weight2 = 1;
		Vec3 L1, L2;
		std::shared_ptr<Light> light = scene.SampleOneLight(&lightSamplingPdf, uLight);
		{
			Vec3 wi;
			Intersection lightPoint;
			real pdf;
			Vec3 Li = light->Sample_Li(isect, &wi, &pdf, &lightPoint, u);
			if (pdf > 0) {
				real scatteringPdf = bsdf->Pdf(isect.wo, wi);
				Vec3 f = bsdf->f(isect.wo, wi);
				VisibilityTester visibilityTester(isect, lightPoint);
				if (!visibilityTester.Unoccluded(scene)) {
					weight1 = PowerHeuristic(1, pdf, 1, scatteringPdf);
					L1 = Li * f * std::abs(isect.n.Dot(wi)) / pdf;
				}

			}

		}
		//Sample BSDF
		{
			Vec3 wi;
			real pdf;
			Vec3 f = bsdf->Sample_f(isect.wo, &wi, &pdf, v);
			Intersection intersection;
			std::shared_ptr<Shape> hitObj;
			if (pdf != 0) {
				real lightPdf = light->Pdf_Li(isect, wi);
				weight2 = PowerHeuristic(1, pdf, 1, lightPdf);
				if (scene.Intersect(isect.SpawnRay(wi), &intersection) && isect.primitive->IsLight()) {
					std::shared_ptr<Light> emissionShape = isect.primitive->GetLight();
					L2 = emissionShape->Emission(intersection, -wi) * f * std::abs(isect.n.Dot(wi)) / pdf;
				}

			}
		}

		return (L1 * weight1 + L2 * weight2) / lightSamplingPdf;
	}

};

NAMESPACE_END