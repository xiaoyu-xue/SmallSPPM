#pragma once

#include "utils.h"
#include "linagl.h"
#include "integrator.h"
#include "sampler_enum.h"
#include "sampler.h"
#include "hashgrid.h"
#include "image_io.h"

#include "debug_utils.h"

NAMESPACE_BEGIN

struct HPoint {
	Vec3 importance, pos, nrm, flux, outDir;
	real r2;
	int64 n; // n = N / ALPHA in the paper
	int64 pix;
	bool used;
	int64 m;
};

class SPPM : public Integrator {
public:
	SPPM(int iterations, int nPhotonsPerStage, int maxDepth, real initialRadius, real alpha, bool batchShrink, const std::shared_ptr<Sampler> &pSampler,
		const std::shared_ptr<SamplerEnum> &pSmplerEnum, bool traceGlossyRay = false) :
		nIterations(iterations), nPhotonsPerRenderStage(nPhotonsPerStage), batchShrink(batchShrink),
		maxDepth(maxDepth), initialRadius(initialRadius), alpha(alpha), sampler(pSampler), samplerEnum(pSmplerEnum), traceGlossyRay(traceGlossyRay)
	{

	}

	//void GeneratePhoton(const Scene &scene, Ray *pr, Vec3 *f, real u, const Vec2 &v, const Vec2 &w) {
	//	real lightPdf;
	//	std::shared_ptr<Light> light = scene.SampleOneLight(&lightPdf, u);
	//	Vec3 Le = light->Emission();
	//	Vec3 pos, lightDir, lightNorm;
	//	real pdfPos, pdfDir;
	//	light->SampleLight(&pos, &lightDir, &lightNorm, &pdfPos, &pdfDir, v, w);
	//	pr->o = pos + lightDir * rayeps;
	//	pr->d = lightDir;
	//	real cosTheta = std::abs(lightNorm.Dot(lightDir));
	//	*f = Le * cosTheta / (pdfPos * pdfDir * lightPdf);
	//}

	void GeneratePhoton(const Scene &scene, Ray *pr, Vec3 *f, real u, const Vec2 &v, const Vec2 &w) {
		real lightPdf;
		std::shared_ptr<Light> light = scene.SampleOneLight(&lightPdf, u);
		//Vec3 Le = light->Emission();
		//Vec3 pos, lightDir, lightNorm;
		Vec3 lightDir;
		Intersection lightPoint;
		real pdfPos, pdfDir;
		Vec3 Le = light->SampleLight(&lightPoint, &lightDir, &pdfPos, &pdfDir, v, w);
		*pr = lightPoint.SpawnRay(lightDir);
		real cosTheta = std::abs(lightPoint.n.Dot(lightDir));
		*f = Le * cosTheta / (pdfPos * pdfDir * lightPdf);
	}

	void TraceEyePath(const Scene &scene, StateSequence &rand, const Ray &ray, int64 pixel) {
		//{
		//	int x = 365, y = 550;
		//	if(pixel == x + scene.GetCamera()->GetFilm()->resX * y) {
		//		debugPixel = 1;
		//		std::cout << "--------------------------------------------------------" << std::endl;
		//	}else {
		//		debugPixel = 0;
		//	}
		//}
		Ray r = ray;
		bool deltaBoundEvent = false;
		Vec3 importance(1.0, 1.0, 1.0);
		for (int i = 0; i < maxDepth; ++i) {
			real t;
			Intersection isect;

			if (!scene.Intersect(r,&isect, &t)) return;
			isect.ComputeScatteringFunction();
			std::shared_ptr<BSDF> bsdf = isect.bsdf;
			Vec3 wi;
			real pdf;

#ifdef _DEBUG
			int primitiveId = isect.primitive->primitiveId;
#endif

			//{
			//	int x = 223, y = 387;
			//	if (pixel == x + scene.GetCamera()->GetFilm()->resX * y) {
			//		//debugPixel = 1;
			//		//std::cout << "isect: " << isect.hit << " " << isect.rayEps << std::endl;
			//		//std::cout << "first hit:" << hitObj->GetId() << std::endl;
			//		//std::cout << isect.hit << "\n" << isect.nl <<"\n" << isect.nl << "\n" << t << std::endl;
			//	}
			//	else {
			//		debugPixel = 0;
			//	}
			//}


			if (bsdf->IsDelta()) {
				Vec3 f = bsdf->Sample_f(-1 * r.d, &wi, &pdf, Vec3(rand(), rand(), rand()));
				importance = f * std::abs(wi.Dot(isect.n)) * importance / pdf;
				//r.o = isect.hit + wi * rayeps + wi.Dot(isect.n) * nEps * isect.n;
				//r.d = wi;
				//r = Ray(isect.hit, wi, Inf, isect.rayEps);
				r = isect.SpawnRay(wi);
				deltaBoundEvent = true;
			}
			else if (traceGlossyRay && bsdf->MatchScatterType(ScatterEventType::BSDF_GLOSSY) && i < maxDepth - 1) {
				//Trace the ray, when bsdf is glossy
				//But this is not standard photon mapping
				HPoint& hp = hitPoints[pixel];
				Vec3 Ld = importance * DirectIllumination(scene, isect, bsdf, rand(), Vec2(rand(), rand()), Vec3(rand(), rand(), rand()));
				directillum[hp.pix] = directillum[hp.pix] + Ld;
				Vec3 f = bsdf->Sample_f(-1 * r.d, &wi, &pdf, Vec3(rand(), rand(), rand()));
				importance = f * std::abs(wi.Dot(isect.n)) * importance / pdf;
				r = isect.SpawnRay(wi);
				deltaBoundEvent = false;
			}
			else {
				HPoint &hp = hitPoints[pixel];
				hp.used = true;
				hp.importance = importance;
				hp.pos = isect.hit;
				hp.nrm = isect.nl;
				hp.pix = pixel;
				hp.outDir = -1 * r.d;
				hitPoints[pixel] = hp;
				if ((i == 0 || deltaBoundEvent) && isect.primitive->IsLight()) {
					std::shared_ptr<Light> emissionShape = isect.primitive->GetLight();
					directillum[hp.pix] = directillum[hp.pix] + importance * emissionShape->Emission(isect, isect.wo);
					//{
					//	if (debugPixel == 1) {
					//		std::cout << "direction illumination(delta): " << hp.importance << " " << 
					//			hitObj->GetEmission() << " " << directillum[hp.pix] << std::endl;
					//	}
					//}
				}
				else {
					Vec3 Ld = hp.importance * DirectIllumination(scene, isect, bsdf, rand(), Vec2(rand(), rand()), Vec3(rand(), rand(), rand()));
					directillum[hp.pix] = directillum[hp.pix] + Ld;
					//{
					//	if (debugPixel == 1) {
					//		std::cout << "direction illumination: " << hp.importance << " " << Ld << " " << directillum[hp.pix] << std::endl;
					//	}
					//}
				}
				return;
			}
		}
	}

	void TracePhoton(const Scene &scene, StateSequence &rand, const Ray &ray, Vec3 photonFlux) {
		Ray r = ray;
		for (int i = 0; i < maxDepth; ++i) {
			real t;
			Intersection isect;
			std::shared_ptr<Shape> hitObj;
			if (!scene.Intersect(r, &isect, &t)) return;

			isect.ComputeScatteringFunction(TransportMode::Importance);
			std::shared_ptr<BSDF> bsdf = isect.bsdf;

			Vec3 wi;
			real pdf;
			Vec3 f = bsdf->Sample_f(-1 * r.d, &wi, &pdf, Vec3(rand(), rand(), rand()));
			Vec3 estimation = f * std::abs(wi.Dot(isect.n)) / pdf;
			if (bsdf->IsDelta()) {
				photonFlux = photonFlux * estimation;
				//r.o = isect.hit + wi * rayeps + wi.Dot(isect.n) * nEps * isect.n;
				//r.d = wi;
				//r = Ray(isect.hit + wi * rayeps + wi.Dot(isect.n) * nEps * isect.n, wi, Inf, isect.rayEps);
				r = isect.SpawnRay(wi);
			}
			else {
				if (i > 0) {
					std::vector<HPoint*> &hp = hashGrid.GetGrid(isect.hit);
					for (HPoint* hitpoint : hp) {
						//Use spinlock, but racing condition is rare when using QMC
						//std::lock_guard<Spinlock> lock(pixelLocks[hitpoint->pix]);
						Vec3 v = hitpoint->pos - isect.hit;
						//if ((hitpoint->nrm.Dot(isect.n) > PhtotonEdgeEps) && (v.Dot(v) <= radius2[hitpoint->pix])) {
						if ((hitpoint->nrm.Dot(isect.nl) > 0.015) && (v.Dot(v) <= radius2[hitpoint->pix])) {
							if (!batchShrink) {
								// unlike N in the paper, hitpoint->n stores "N / ALPHA" to make it an integer value
								real g = (photonNums[hitpoint->pix] * alpha + alpha) / (photonNums[hitpoint->pix] * alpha + 1.f);
								radius2[hitpoint->pix] = radius2[hitpoint->pix] * g;
								photonNums[hitpoint->pix] += 1;
								Vec3 contribution = hitpoint->importance * bsdf->f(hitpoint->outDir, -1 * r.d) * photonFlux;
								flux[hitpoint->pix] = (flux[hitpoint->pix] + contribution) * g;
							}
							else {
								hitpoint->m++;
								Vec3 contribution = hitpoint->importance * bsdf->f(hitpoint->outDir, -1 * r.d) * photonFlux;
								flux[hitpoint->pix] = (flux[hitpoint->pix] + contribution);
							}
						}
					}
				}
				//real p = estimation.maxValue();
				real p = std::min((real)1.0, (estimation * photonFlux).Y() / photonFlux.Y());
				if (p < 1) {
					if (rand() < p) {
						photonFlux = photonFlux / p;
					}
					else {
						break;
					}
				}
				photonFlux = photonFlux * estimation;
				//r.o = isect.hit + wi * rayeps + wi.Dot(isect.n) * nEps * isect.n;
				//r.d = wi;
				//r = Ray(isect.hit + wi * rayeps + wi.Dot(isect.n) * nEps * isect.n, wi, Inf, isect.rayEps);
				r = isect.SpawnRay(wi);
			}
		}
	}

	void GenerateRadiusImage(const Scene &scene) {
		int resX = scene.GetCamera()->GetFilm()->resX;
		int resY = scene.GetCamera()->GetFilm()->resY;
		std::vector<Vec3> radImg(resX * resY);
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
				real val = 1.f - (std::sqrt(radius2[index]) - std::sqrt(minRadius2)) /
					(std::sqrt(maxRadius2) - std::sqrt(minRadius2));
				radImg[index].x = val;
				radImg[index].y = val;
				radImg[index].z = val;
			}
		}

		ImageIO::WriteBmpFile("rad_image.bmp", radImg, resX, resY, 2.2f);
	}

	void Render(const Scene &scene) override {
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
					Vec2 pixelSample(u, v);
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
					hashGrid.AddPoint(std::move(std::pair<Vec3, HPoint*>(hp.pos, &hp)), std::sqrt(radius2[i]));
				}
			}
			hashGrid.BuildHashGrid(std::sqrt(maxRadius2) + eps);

			//Trace photon
			ParallelFor((int64)0, nPhotonsPerRenderStage, [&](int64 j) {
				Ray ray;
				Vec3 photonFlux;
				RandomStateSequence rand(sampler, iter * nPhotonsPerRenderStage + j);
				GeneratePhoton(scene, &ray, &photonFlux, rand(), Vec2(rand(), rand()), Vec2(rand(), rand()));
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

			real percentage = 100.f * (iter + 1) / nIterations;
			fprintf(stderr, "\rIterations: %5.2f%%", percentage);
		}

		// density estimation
		std::cout << "\nFlux size: " << flux.size() << std::endl;
		ParallelFor(0, (int)flux.size(), [&](int i) {
			c[i] = c[i] + flux[i] * (1.f / (PI * radius2[i] * nIterations * nPhotonsPerRenderStage))
				+ directillum[i] / (real)nIterations;
			//if (std::isnan(c[i].x)  || std::isnan(c[i].y)  || std::isnan(c[i].z)) {
			//	std::cout << i % resX << " " << i / resX << " "<< i << std::endl;
			//}
			//c[i] = c[i] + directillum[i] / nIterations;
		});
		//{
		//	for (int x = 356; x <= 374; ++x) {
		//		for (int y = 545; y <= 560; ++y) {
		//			int index = x + y * resX;
		//			if (c[index].x >= 1 || c[index].y >= 1 || c[index].z >= 1) {
		//				std::cout << c[index] << std::endl;
		//				std::cout << x << " " << y << std::endl;
		//			}

		//				
		//		}
		//	}
		//	int x = 365, y = 550;
		//	int index = x + y * resX;
		//	std::cout << c[index] << std::endl;
		//}
		//std::cout << flux[591714] << " " << radius2[591714] << std::endl;
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
			flux[i] = Vec3();
			hitPoints[i].m = 0;
		}
	}

	HashGrid<HPoint*> hashGrid;
	std::shared_ptr<Sampler> sampler;
	std::shared_ptr<SamplerEnum> samplerEnum;
	std::vector<real> radius2;
	//std::vector<int64> photonNums;
	std::vector<real> photonNums;
	std::vector<Vec3> flux;
	std::vector<Vec3> directillum;
	std::vector<HPoint> hitPoints;
	std::vector<Vec3> c;
	std::vector<Spinlock> pixelLocks;
	const real initialRadius;
	const int maxDepth;
	const int nIterations;
	const int64 nPhotonsPerRenderStage;
	const real alpha;
	bool batchShrink;
	bool traceGlossyRay;
};

NAMESPACE_END