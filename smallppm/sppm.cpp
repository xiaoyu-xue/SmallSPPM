#include "sppm.h"
#include "image_io.h"
#include "scene.h"
#include "threading.h"
NAMESPACE_BEGIN

void SPPM::GeneratePhoton(const Scene& scene, Ray* pr, Vec3* f, real u, const Vec2& v, const Vec2& w) {
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
	if (pdfPos == 0 || pdfDir == 0) {
		*f = 0;
	}
	else {
		*f = Le * cosTheta / (pdfPos * pdfDir * lightPdf);
	}
}

void SPPM::TraceEyePath(const Scene& scene, StateSequence& rand, const Ray& ray, int64 pixel, MemoryArena& arena) {
	Ray r = ray;
	bool deltaBoundEvent = false;
	Vec3 importance(1.0, 1.0, 1.0);
	for (int i = 0; i < maxDepth; ++i) {
		Intersection isect;
		if (!scene.Intersect(r, &isect)) {
			//Environment 
			if ((i == 0 || deltaBoundEvent) && scene.GetEnvironmentLight()) {
				directillum[pixel] += importance * scene.GetEnvironmentLight()->Emission(r);
			}
			return;
		}
		scene.QueryIntersectionInfo(r, &isect);
		isect.ComputeScatteringFunction(arena);
		BSDF* bsdf = isect.bsdf;
		Vec3 wi;
		real pdf;

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
			HPoint& hp = hitPoints[pixel];
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

			}
			else {
				Vec3 Ld = hp.importance * DirectIllumination(scene, isect, bsdf, rand(), Vec2(rand(), rand()), Vec3(rand(), rand(), rand()));
				directillum[hp.pix] = directillum[hp.pix] + Ld;

			}
			return;
		}
	}
}

void SPPM::TracePhoton(const Scene& scene, StateSequence& rand, const Ray& ray, Vec3 photonFlux, MemoryArena& arena) {
	Ray r = ray;
	for (int i = 0; i < maxDepth; ++i) {
		Intersection isect;
		if (!scene.Intersect(r, &isect)) return;
		scene.QueryIntersectionInfo(r, &isect);
		isect.ComputeScatteringFunction(arena, TransportMode::Importance);
		BSDF* bsdf = isect.bsdf;

		Vec3 wi;
		real pdf;
		Vec3 f = bsdf->Sample_f(-1 * r.d, &wi, &pdf, Vec3(rand(), rand(), rand()));
		Vec3 estimation = f * std::abs(wi.Dot(isect.n)) / pdf;
		if (bsdf->IsDelta()) {
			//photonFlux = photonFlux * estimation;
			//r = isect.SpawnRay(wi);
		}
		else {
			if (i > 0) {
				std::vector<HPoint*>& hp = hashGrid.GetGrid(isect.hit);
				for (HPoint* hitpoint : hp) {
					//Use spinlock, but racing condition is rare when using QMC
					//std::lock_guard<Spinlock> lock(pixelLocks[hitpoint->pix]);
					Vec3 v = hitpoint->pos - isect.hit;
					//if ((hitpoint->nrm.Dot(isect.n) > PhtotonEdgeEps) && (v.Dot(v) <= radius2[hitpoint->pix])) {
					if ((hitpoint->nrm.Dot(isect.nl) > 0.0015) && (v.Dot(v) <= radius2[hitpoint->pix])) {
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
		r = isect.SpawnRay(wi);
	}
}

void SPPM::GenerateRadiusImage(const Scene& scene) {
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

void SPPM::Render(const Scene& scene) {
	fprintf(stderr, "Rendering ...\n");
	int resX = scene.GetCamera()->GetFilm()->resX;
	int resY = scene.GetCamera()->GetFilm()->resY;
	Initialize(resX, resY);
	for (int iter = 0; iter < nIterations; ++iter) {
		//Trace eye_path
		std::vector<MemoryArena> memoryArenas(ThreadsNumber());
		ParallelFor(0, resY, [&](int y) {
			for (int x = 0; x < resX; x++) {
				MemoryArena& arena = memoryArenas[ThreadIndex()];
				int pixel = x + y * resX;
				hitPoints[pixel].used = false;
				uint64 instance = samplerEnum->GetIndex(iter, x, y);
				RandomStateSequence rand(sampler, instance);
				real u = samplerEnum->SampleX(x, rand());
				real v = samplerEnum->SampleY(y, rand());
				Vec2 pixelSample(u, v);
				Ray ray = scene.GetCamera()->GenerateRay(x, y, pixelSample);
				TraceEyePath(scene, rand, ray, pixel, arena);
			}
		});

		hashGrid.ClearHashGrid();
		real maxRadius2 = 0.0;
		for (int i = 0; i < (int)hitPoints.size(); ++i) {
			HPoint& hp = hitPoints[i];
			if (hp.used) {
				maxRadius2 = std::max(maxRadius2, radius2[i]);
				hashGrid.AddPoint(std::move(std::pair<Vec3, HPoint*>(hp.pos, &hp)), std::sqrt(radius2[i]));
			}
		}
		hashGrid.BuildHashGrid(std::sqrt(maxRadius2) + eps);

		//Trace photon
		ParallelFor((int64)0, nPhotonsPerRenderStage, [&](int64 j) {
			MemoryArena& arena = memoryArenas[ThreadIndex()];
			Ray ray;
			Vec3 photonFlux;
			RandomStateSequence rand(sampler, iter * nPhotonsPerRenderStage + j);
			GeneratePhoton(scene, &ray, &photonFlux, rand(), Vec2(rand(), rand()), Vec2(rand(), rand()));
			TracePhoton(scene, rand, ray, photonFlux, arena);
			arena.Reset();
		});

		//Update flux, radius, photonNums if batchShrink
		if (batchShrink) {
			size_t nHitPoints = hitPoints.size();
			ParallelFor(size_t(0), nHitPoints, [&](size_t i) {
				HPoint& hp = hitPoints[i];
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
		//c[i] = c[i] + flux[i] * (1.f / (PI * radius2[i] * nIterations * nPhotonsPerRenderStage));
	});

	scene.GetCamera()->GetFilm()->SetImage(c);
	//GenerateRadiusImage(scene);
}

NAMESPACE_END