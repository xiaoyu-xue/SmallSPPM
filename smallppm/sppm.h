#pragma once

#include "utils.h"
#include "linagl.h"
#include "integrator.h"
#include "sampler_enum.h"
#include "sampler.h"
#include "hashgrid.h"
#include "image_io.h"

NAMESPACE_BEGIN

struct HPoint {
	Vec importance, pos, nrm, flux, outDir;
	real r2;
	int64 n; // n = N / ALPHA in the paper
	int64 pix;
	bool used;
	int64 m;
};

class SPPM : public Integrator {
public:
	SPPM(int iterations, int nPhotonsPerStage, int maxDepth, real initialRadius, real alpha, bool batchShrink, const std::shared_ptr<Sampler> &pSampler,
		const std::shared_ptr<SamplerEnum> &pSmplerEnum) :
		nIterations(iterations), nPhotonsPerRenderStage(nPhotonsPerStage), batchShrink(batchShrink),
		maxDepth(maxDepth), initialRadius(initialRadius), alpha(alpha), sampler(pSampler), samplerEnum(pSmplerEnum)
	{

	}

	void GeneratePhoton(const Scene &scene, Ray *pr, Vec *f, real u, const Vec &v, const Vec &w) {
		real lightPdf;
		std::shared_ptr<Light> light = scene.SampleOneLight(&lightPdf, u);
		Vec Le = light->Emission();
		Vec pos, lightDir, lightNorm;
		real pdfPos, pdfDir;
		light->SampleLight(&pos, &lightDir, &lightNorm, &pdfPos, &pdfDir, v, w);
		pr->o = pos + lightDir * rayeps;
		pr->d = lightDir;
		real cosTheta = std::abs(lightNorm.dot(lightDir));
		*f = Le * cosTheta / (pdfPos * pdfDir * lightPdf);
	}

	void TraceEyePath(const Scene &scene, StateSequence &rand, const Ray &ray, int64 pixel) {
		Ray r = ray;
		bool deltaBoundEvent = false;
		Vec importance(1.0, 1.0, 1.0);
		for (int i = 0; i < maxDepth; ++i) {
			real t;
			Intersection isect;
			std::shared_ptr<Shape> hitObj;
			if (!scene.Intersect(r, &t, &isect, hitObj)) return;
			std::shared_ptr<BSDF> bsdf = hitObj->GetBSDF(isect);
			Vec wi;
			real pdf;
			if (bsdf->IsDelta()) {
				Vec f = bsdf->Sample_f(-1 * r.d, &wi, &pdf, Vec(rand(), rand(), rand()));
				importance = f * std::abs(wi.dot(isect.n)) * importance / pdf;
				r.o = isect.hit + wi * rayeps + wi.dot(isect.n) * nEps * isect.n;
				r.d = wi;
				deltaBoundEvent = true;
			}
			else {
				HPoint &hp = hitPoints[pixel];
				hp.used = true;
				hp.importance = importance;
				hp.pos = isect.hit;
				hp.nrm = isect.n;
				hp.pix = pixel;
				hp.outDir = -1 * r.d;
				hitPoints[pixel] = hp;
				if ((i == 0 || deltaBoundEvent) && hitObj->IsLight()) {
					directillum[hp.pix] = directillum[hp.pix] + importance * hitObj->GetEmission();
				}
				else {
					Vec Ld = hp.importance * DirectIllumination(scene, isect, bsdf, rand(), Vec(rand(), rand(), rand()), Vec(rand(), rand(), rand()));
					directillum[hp.pix] = directillum[hp.pix] + Ld;
				}
				return;
			}
		}
	}

	void TracePhoton(const Scene &scene, StateSequence &rand, const Ray &ray, Vec photonFlux) {
		Ray r = ray;
		for (int i = 0; i < maxDepth; ++i) {
			real t;
			Intersection isect;
			std::shared_ptr<Shape> hitObj;
			if (!scene.Intersect(r, &t, &isect, hitObj)) return;
			std::shared_ptr<BSDF> bsdf = hitObj->GetBSDF(isect, TransportMode::Importance);
			Vec wi;
			real pdf;
			Vec f = bsdf->Sample_f(-1 * r.d, &wi, &pdf, Vec(rand(), rand(), rand()));
			Vec estimation = f * std::abs(wi.dot(isect.n)) / pdf;
			if (bsdf->IsDelta()) {
				photonFlux = photonFlux * estimation;
				r.o = isect.hit + wi * rayeps + wi.dot(isect.n) * nEps * isect.n;
				r.d = wi;
			}
			else {
				if (i > 0) {
					std::vector<HPoint*> &hp = hashGrid.GetGrid(isect.hit);
					for (HPoint* hitpoint : hp) {
						//Use spinlock, but racing condition is rare when using QMC
						//std::lock_guard<Spinlock> lock(pixelLocks[hitpoint->pix]);
						Vec v = hitpoint->pos - isect.hit;
						if ((hitpoint->nrm.dot(isect.n) > PhtotonEdgeEps) && (v.dot(v) <= radius2[hitpoint->pix])) {
							if (!batchShrink) {
								// unlike N in the paper, hitpoint->n stores "N / ALPHA" to make it an integer value
								real g = (photonNums[hitpoint->pix] * alpha + alpha) / (photonNums[hitpoint->pix] * alpha + 1.0);
								radius2[hitpoint->pix] = radius2[hitpoint->pix] * g;
								photonNums[hitpoint->pix] += 1;
								Vec contribution = hitpoint->importance * bsdf->f(hitpoint->outDir, -1 * r.d) * photonFlux;
								flux[hitpoint->pix] = (flux[hitpoint->pix] + contribution) * g;
							}
							else {
								hitpoint->m++;
								Vec contribution = hitpoint->importance * bsdf->f(hitpoint->outDir, -1 * r.d) * photonFlux;
								flux[hitpoint->pix] = (flux[hitpoint->pix] + contribution);
							}
						}
					}
				}
				//real p = estimation.maxValue();
				real p = std::min(1.0, (estimation * photonFlux).Y() / photonFlux.Y());
				if (p < 1) {
					if (rand() < p) {
						photonFlux = photonFlux / p;
					}
					else {
						break;
					}
				}
				photonFlux = photonFlux * estimation;
				r.o = isect.hit + wi * rayeps + wi.dot(isect.n) * nEps * isect.n;
				r.d = wi;
			}
		}
	}

	void GenerateRadiusImage(const Scene &scene) {
		int resX = scene.GetCamera()->GetFilm()->resX;
		int resY = scene.GetCamera()->GetFilm()->resY;
		std::vector<Vec> radImg(resX * resY);
		real minRadius2 = Inf, maxRadius2 = 0.0, avgRadius = 0.0;
		for (int y = 0; y < resY; ++y) {
			for (int x = 0; x < resX; ++x) {
				int index = x + y * resX;
				minRadius2 = std::min(minRadius2, radius2[index]);
				maxRadius2 = std::max(maxRadius2, radius2[index]);
				avgRadius += std::sqrt(radius2[index]);
			}
		}
		std::cout << "minRasius: " << std::sqrt(minRadius2) << std::endl;
		std::cout << "avgRasius: " << avgRadius / resX / resY << std::endl;
		for (int y = 0; y < resY; ++y) {
			for (int x = 0; x < resX; ++x) {
				int index = x + y * resX;
				real val = 1.0 - (std::sqrt(radius2[index]) - std::sqrt(minRadius2)) /
					(std::sqrt(maxRadius2) - std::sqrt(minRadius2));
				radImg[index].x = val;
				radImg[index].y = val;
				radImg[index].z = val;
			}
		}

		ImageIO::WriteBmpFile("rad_image.bmp", radImg, resX, resY, 2.2);
	}

	void Render(const Scene &scene) {
		fprintf(stderr, "Rendering ...\n");
		int resX = scene.GetCamera()->GetFilm()->resX;
		int resY = scene.GetCamera()->GetFilm()->resY;
		Initialize(resX, resY);
		for (int iter = 0; iter < nIterations; ++iter) {
			//Trace eye_path
			ParallelFor(0, resY, [&](int y) {
				for (int x = 0; x < resX; x++) {
					int pixel = x + y * resX;
					hitPoints[pixel].used = false;
					uint64 instance = samplerEnum->GetIndex(iter, x, y);
					RandomStateSequence rand(sampler, instance);
					real u = samplerEnum->SampleX(x, rand());
					real v = samplerEnum->SampleY(y, rand());
					Vec pixelSample(u, v, 0);
					Ray ray = scene.GetCamera()->GenerateRay(x, y, pixelSample);
					TraceEyePath(scene, rand, ray, pixel);
				}
			});

			hashGrid.ClearHashGrid();
			real maxRadius2 = 0.0;
			for (int i = 0; i < (int)hitPoints.size(); ++i) {
				HPoint &hp = hitPoints[i];
				if (hp.used) {
					maxRadius2 = std::max(maxRadius2, radius2[i]);
					hashGrid.AddPoint(std::move(std::pair<Vec, HPoint*>(hp.pos, &hp)), std::sqrt(radius2[i]));
				}
			}
			hashGrid.BuildHashGrid(std::sqrt(maxRadius2) + eps);

			//Trace photon
			ParallelFor((int64)0, nPhotonsPerRenderStage, [&](int64 j) {
				Ray ray;
				Vec photonFlux;
				RandomStateSequence rand(sampler, iter * nPhotonsPerRenderStage + j);
				GeneratePhoton(scene, &ray, &photonFlux, rand(), Vec(rand(), rand(), rand()), Vec(rand(), rand(), rand()));
				TracePhoton(scene, rand, ray, photonFlux);
			});

			//Update flux, radius, photonNums if batchShrink
			if(batchShrink) {
				size_t nHitPoints = hitPoints.size();
				ParallelFor(size_t(0), nHitPoints, [&](size_t i) {
					HPoint &hp = hitPoints[i];
					if (hp.m > 0) {
						real g = (photonNums[hp.pix] + alpha * hp.m) / (photonNums[hp.pix] + hp.m);
						radius2[hp.pix] = radius2[hp.pix] * g;
						flux[hp.pix] = flux[hp.pix] * g;
						photonNums[hp.pix] = photonNums[hp.pix] + alpha * hp.m;
						hp.m = 0;
					}
				});
			}

			real percentage = 100.*(iter + 1) / nIterations;
			fprintf(stderr, "\rIterations: %5.2f%%", percentage);
		}

		// density estimation
		std::cout << "\nFlux size: " << flux.size() << std::endl;
		ParallelFor(0, (int)flux.size(), [&](int i) {
			c[i] = c[i] + flux[i] * (1.0 / (PI * radius2[i] * nIterations * nPhotonsPerRenderStage))
				+ directillum[i] / nIterations;
			//c[i] = c[i] + directillum[i] / nIterations;
		});
		scene.GetCamera()->GetFilm()->SetImage(c);
		//GenerateRadiusImage(scene);
	}

private:
	void Initialize(int w, int h) {
		radius2.resize(w * h);
		photonNums.resize(w * h);
		flux.resize(w * h);
		hitPoints.resize(w * h);
		directillum.resize(w * h);
		c.resize(w * h);
		pixelLocks.resize(w * h);
		for (int i = 0; i < w * h; ++i) {
			radius2[i] = initialRadius * initialRadius;
			photonNums[i] = 0;
			flux[i] = Vec();
			hitPoints[i].m = 0;
		}
	}

	HashGrid<HPoint*> hashGrid;
	std::shared_ptr<Sampler> sampler;
	std::shared_ptr<SamplerEnum> samplerEnum;
	std::vector<real> radius2;
	//std::vector<int64> photonNums;
	std::vector<real> photonNums;
	std::vector<Vec> flux;
	std::vector<Vec> directillum;
	std::vector<HPoint> hitPoints;
	std::vector<Vec> c;
	std::vector<Spinlock> pixelLocks;
	const real initialRadius;
	const int maxDepth;
	const int nIterations;
	const int64 nPhotonsPerRenderStage;
	const real alpha;
	bool batchShrink;
};

NAMESPACE_END