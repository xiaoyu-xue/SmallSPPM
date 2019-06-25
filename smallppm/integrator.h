#pragma once

#include "utils.h"
#include "scene.h"
#include "intersection.h"
#include "visibility.h"


real PowerHeuristic(int nf, real fPdf, int ng, real gPdf) {
	real f = nf * fPdf, g = ng * gPdf;
	return (f * f) / (f * f + g * g);
}

class Integrator {
public:
	virtual void Render(const Scene &scene) = 0;
	virtual ~Integrator() {}
	static Vec DirectIllumination(const Scene &scene, const Intersection &isect, const std::shared_ptr<BSDF> &bsdf,
		real uLight, Vec u, Vec v) {
		Vec L;
		/*
		const std::vector<std::shared_ptr<Light>> lights = scene.GetLights();
		for (auto light : lights) {
			Vec dir;
			std::shared_ptr<Shape> hitObj;
			real t;
			Intersection lightPoint;
			Vec Li = light->DirectIllumination(isect, bsdf, importance, &dir, &lightPoint, v);
			Intersection intersection;
			if (scene.Intersect(Ray(isect.hit, dir), &t, &intersection, hitObj) && hitObj->GetId() == light->GetId()) {
				L = L + Li;
			}
		}
		return L;
		*/

		//Sample light
		real lightSamplingPdf;
		real weight1 = 1, weight2 = 1;
		Vec L1, L2;
		std::shared_ptr<Light> light = scene.SampleOneLight(&lightSamplingPdf, uLight);
		{
			Vec wi;
			Intersection lightPoint;
			real pdf;
			Vec Li = light->Sample_Li(isect, &wi, &pdf, &lightPoint, u);
			if (pdf > 0) {
				real scatteringPdf = bsdf->Pdf(isect.wo, wi);
				Vec f = bsdf->f(isect.wo, wi);
				VisibilityTester visibilityTester(isect, lightPoint);
				if (!visibilityTester.Unoccluded(scene)) {
					weight1 = PowerHeuristic(1, pdf, 1, scatteringPdf);
					L1 = Li * f * std::abs(isect.n.dot(wi)) / pdf;
				}
			}

		}


		//Sample BSDF
		{
			Vec wi;
			real pdf;
			Vec f = bsdf->Sample_f(isect.wo, &wi, &pdf, v);
			Intersection intersection;
			std::shared_ptr<Shape> hitObj;
			real t;
			if (pdf != 0) {
				real lightPdf = light->Pdf_Li(isect, wi);
				weight2 = PowerHeuristic(1, pdf, 1, lightPdf);
				if (scene.Intersect(Ray(isect.hit + wi * rayeps, wi), &t, &intersection, hitObj) && hitObj->IsLight()) {
					L2 = hitObj->GetEmission() * f * std::abs(isect.n.dot(wi)) / pdf;
				}
			}
		}
		//return L1 / lightSamplingPdf;
		//return L2 / lightSamplingPdf;
		return  (L1 * weight1 + L2 * weight2) / lightSamplingPdf;

	}
};