#pragma once

#include "BidirectionalRenderer.h"
#include "core/Scene.h"
#include "core/Intersection.h"
#include "light/Light.h"

GYT_NAMESPACE_BEGIN

int BidirectionalRenderer::GenerateLightPath(
	const Scene				&scene, 
	StateSequence			&rand,
	MemoryPool				&arena, 
	std::vector<PathVertex>	&lightPath, 
	int						maxDepth)
{
	real pdfLight;
	real pdfA;
	real pdfDir;
	Vec3 Le;
	Vec3 dir;
	Intersection isect;
	Light* light = scene.SampleOneLight(&pdfLight, rand());
	if (light->IsAreaLight()) {
		Le = light->Emission();
		light->SampleOnLight(&isect, &dir, &pdfA, &pdfDir, Vec2(rand(), rand()), Vec2(rand(), rand()));
	}
	else {
		Le = light->Emission();
		light->SampleLight(&isect, &dir, &pdfA, &pdfDir, Vec2(rand(), rand()), Vec2(rand(), rand()));
	}
	lightPath[0].mIsect = isect;
	lightPath[0].mPdfFwd = pdfA;
	lightPath[0].mThroughput = Le;
	lightPath[0].mIsDelta = light->IsDeltaLight();
	real cosTheta = std::abs(Dot(isect.mNormal, dir));
	Vec3 throughput = Le * cosTheta / pdfLight / pdfA / pdfDir;
	Ray ray(isect.mPos + dir * RayEps, dir);

	int nVertices = BidirectionalRenderer::Trace(scene, arena, rand, ray, 1, throughput, pdfDir, lightPath, maxDepth, TransportMode::Importance);
	return nVertices;
}

int BidirectionalRenderer::GenerateCameraPath(
	const Scene				&scene, 
	const Camera			&camera, 
	StateSequence			&rand, 
	MemoryPool				&arena, 
	std::vector<PathVertex>	&cameraPath, 
	const Ray				&cameraRay, 
	int						maxdDpth)
{
	if (maxdDpth == 0) return 0;
	Vec3 throughput(1.f, 1.f, 1.f);
	cameraPath[0].mIsect.mPos = camera.mPos;
	cameraPath[0].mIsect.mNormal = camera.mCz;
	cameraPath[0].mPdfFwd = camera.PdfPos();
	cameraPath[0].mThroughput = throughput;
	real pdfDir = camera.PdfDir(cameraRay);
	const Camera* pCam = &camera;
	cameraPath[0].mpCamera = const_cast<Camera*>(pCam);
	int nCameraVertices = BidirectionalRenderer::Trace(scene, arena, rand, cameraRay, 1, throughput, pdfDir, cameraPath, maxdDpth, TransportMode::Radiance);
	return nCameraVertices;
}

int BidirectionalRenderer::Trace(
	const Scene				&scene,
	MemoryPool				&arena,
	StateSequence			&rand,
	const Ray				&r,
	int						depth,
	Vec3					throughput,
	real					pdfFwd,
	std::vector<PathVertex>	&path,
	int						maxDepth,
	TransportMode			mode)
{

	int bounce = 1;
	Ray ray = r;
	real pdfW = pdfFwd;
	while (true) {
		PathVertex& prev = path[bounce - 1];
		PathVertex& vertex = path[bounce];
		Intersection& isect = path[bounce].mIsect;
		if (!scene.Intersect(ray, &isect)) {
			break;
		}
		scene.QueryIntersectionInfo(ray, &isect);
		isect.ComputeScatteringFunction(arena);
		BSDF* pBSDF = isect.mpBSDF;

		path[bounce].mPdfFwd = ConvertSolidToArea(pdfW, prev, vertex);
		path[bounce].mThroughput = throughput;
		path[bounce].mIsDelta = pBSDF->IsDelta();

		Vec3 wi;
		Vec3 f = pBSDF->Sample(-ray.mDir, &wi, &pdfW, Vec3(rand(), rand(), rand()));
		Vec3 estimation;

		estimation = f * std::abs(Dot(isect.mNormal, wi)) / pdfW;
		throughput = throughput * estimation;

		ray = isect.SpawnRay(wi);

		real pdfWPrev = isect.mpBSDF->Pdf(wi, -1 * ray.mDir);
		prev.mPdfPrev = ConvertSolidToArea(pdfWPrev, vertex, prev);

		bounce++;
		if (bounce >= maxDepth + 1) break;
	}
	return bounce;
}

real BidirectionalRenderer::ConvertSolidToArea(
	real				pdfW, 
	const PathVertex	&vertex, 
	const PathVertex	&nxt)
{
	Vec3 dir = nxt.mIsect.mPos - vertex.mIsect.mPos;
	real dist2 = dir.Length2();
	if (dist2 == 0) return 0;
	dir.Normalize();
	real cosTheta = std::abs(dir.Dot(nxt.mIsect.mNormal));
	return pdfW * cosTheta / dist2;
}

bool BidirectionalRenderer::IsConnectable(
	const Scene			&scene, 
	const Vec3			&pointA, 
	const Vec3			&pointB)
{
	Vec3 dir = pointA - pointB;
	real dist = dir.Length();
	Ray ray(pointB, dir.Norm(), 0.0, dist - Eps);  //Be careful ! Float error will cause the incorrect intersection
	if (scene.Intersect(ray)) return false;
	return true;
}

real BidirectionalRenderer::G(
	const PathVertex		&vertexA, 
	const PathVertex		&vertexB)
{
	Vec3 dirAtoB = vertexB.mIsect.mPos - vertexA.mIsect.mPos;
	real dist2 = dirAtoB.Length2();
	dirAtoB.Norm();
	real cosThetaA = std::abs(dirAtoB.Dot(vertexA.mIsect.mNormal));
	real cosThetaB = std::abs(dirAtoB.Dot(vertexB.mIsect.mNormal));
	real geometryTerm = cosThetaA * cosThetaB / dist2;
	return geometryTerm;
}

real BidirectionalRenderer::MISWeight(
	const Scene				&scene,
	StateSequence			&rand,
	std::vector<PathVertex>	&lightPath, 
	std::vector<PathVertex>	&cameraPath, 
	int						s, 
	int						t, 
	PathVertex				&sampled)
{
	if (s + t == 2) return 1;

	PathVertex* lightVertex = s > 0 ? &lightPath[s - 1] : nullptr,
		* cameraVertex = t > 0 ? &cameraPath[t - 1] : nullptr,
		* lightVertexMinus = s > 1 ? &lightPath[s - 2] : nullptr,
		* cameraVertexMinus = t > 1 ? &cameraPath[t - 2] : nullptr;

	ScopedAssignment<PathVertex> a1;
	if (s == 1)
		a1 = { lightVertex, sampled };
	else if (t == 1)
		a1 = { cameraVertex, sampled };


	ScopedAssignment<bool> a2, a3;
	if (lightVertex) a2 = { &lightVertex->mIsect.mIsDelta, false };
	if (cameraVertex) a3 = { &cameraVertex->mIsect.mIsDelta, false };

	ScopedAssignment<real> a4;
	if (cameraVertex) {
		if (s > 0) {
			real pdfW, pdfA;
			if (lightVertexMinus == nullptr) {
				Vec3 lightToHitPoint = cameraVertex->mIsect.mPos - lightVertex->mIsect.mPos;
				lightToHitPoint.Normalize();
				real cosTheta = lightVertex->mIsect.mNormal.Dot(lightToHitPoint);
				pdfW = CosineHemispherePdf(cosTheta);
				pdfA = ConvertSolidToArea(pdfW, *lightVertex, *cameraVertex);
			}
			else {
				Vec3 wo = (lightVertexMinus->mIsect.mPos - lightVertex->mIsect.mPos).Norm();
				Vec3 wi = (cameraVertex->mIsect.mPos - lightVertex->mIsect.mPos).Norm();
				pdfW = lightVertex->mIsect.mpBSDF->Pdf(wo, wi);
				pdfA = ConvertSolidToArea(pdfW, *lightVertex, *cameraVertex);
			}
			a4 = { &cameraVertex->mPdfPrev, pdfA };
		}
		else {

			real pdfLight, pdfPos;
			const Light* pLight = scene.SampleOneLight(&pdfLight, rand());
			Intersection isect;
			pLight->SampleOnePoint(&isect, &pdfPos, Vec2(rand(), rand()));
			a4 = { &cameraVertex->mPdfPrev, pdfLight * pdfPos };
		}
	}

	ScopedAssignment<real> a5;
	if (cameraVertexMinus) {
		if (s > 0) {
			Vec3 wo = (lightVertex->mIsect.mPos - cameraVertex->mIsect.mPos).Norm();
			Vec3 wi = (cameraVertexMinus->mIsect.mPos - cameraVertex->mIsect.mPos).Norm();
			real pdfW = cameraVertex->mIsect.mpBSDF->Pdf(wo, wi);
			real pdfA = ConvertSolidToArea(pdfW, *cameraVertex, *cameraVertexMinus);
			a5 = { &cameraVertexMinus->mPdfPrev, pdfA };
		}
		else {
			real pdfLight;
			Light* pLight = scene.SampleOneLight(&pdfLight, rand());
			Vec3 hitPointToLight = cameraVertex->mIsect.mPos - cameraVertexMinus->mIsect.mPos;
			real dist = hitPointToLight.Length();
			hitPointToLight.Normalize();
			Vec3 lightToHitPoint = -1 * hitPointToLight;
			Vec3 lightNormal = cameraVertex->mIsect.mNormal;
			real cosThetaLight = std::abs(lightNormal.Dot(lightToHitPoint));
			real pdfW = CosineHemispherePdf(cosThetaLight);
			real pdfA = pdfW * (std::abs(cameraVertexMinus->mIsect.mNormal.Dot(hitPointToLight))) / dist / dist;
			a5 = { &cameraVertexMinus->mPdfPrev, pdfA * pdfLight};
		}

	}

	// Update reverse density of vertices $\pq{}_{s-1}$ and $\pq{}_{s-2}$
	ScopedAssignment<real> a6;
	if (lightVertex) {
		real pdfW, pdfA;
		if (cameraVertexMinus == nullptr) {
			Vec3 cameraToHitPoint = lightVertex->mIsect.mPos - cameraVertex->mIsect.mPos;
			cameraToHitPoint.Normalize();
			real cosTheata = std::abs(cameraVertex->mIsect.mNormal.Dot(cameraToHitPoint));
			Ray ray(cameraVertex->mIsect.mPos, cameraToHitPoint);
			pdfW = cameraVertex->mpCamera->PdfDir(ray);
			pdfA = ConvertSolidToArea(pdfW, *cameraVertex, *lightVertex);
		}
		else {
			Vec3 wo = (cameraVertexMinus->mIsect.mPos - cameraVertex->mIsect.mPos).Norm();
			Vec3 wi = (lightVertex->mIsect.mPos - cameraVertex->mIsect.mPos).Norm();
			pdfW = cameraVertex->mIsect.mpBSDF->Pdf(wo, wi);
			pdfA = ConvertSolidToArea(pdfW, *cameraVertex, *lightVertex);
		}
		a6 = { &lightVertex->mPdfPrev, pdfA };
	}

	ScopedAssignment<real> a7;
	if (lightVertexMinus) {
		Vec3 wo = (cameraVertex->mIsect.mPos - lightVertex->mIsect.mPos).Norm();
		Vec3 wi = (lightVertexMinus->mIsect.mPos - lightVertex->mIsect.mPos).Norm();
		real PdfW = lightVertex->mIsect.mpBSDF->Pdf(wo, wi);
		real PdfA = ConvertSolidToArea(PdfW, *lightVertex, *lightVertexMinus);
		a7 = { &lightVertexMinus->mPdfPrev, PdfA };
	}


	real sumRi = 0.0;
	auto remap0 = [](real f)->real {return f != 0 ? f : 1; };

	std::vector<real> cc;

	real ri = 1.0;
	for (int i = t - 1; i > 0; --i) {
		ri *= remap0(cameraPath[i].mPdfPrev) / remap0(cameraPath[i].mPdfFwd);
		if (!cameraPath[i].mIsect.mIsDelta && !cameraPath[i - 1].mIsect.mIsDelta)
			sumRi += ri;
	}

	std::vector<real> dd;
	ri = 1.0;
	for (int i = s - 1; i >= 0; --i) {
		ri *= remap0(lightPath[i].mPdfPrev) / remap0(lightPath[i].mPdfFwd);
		bool deltalightVertex = i > 0 ? lightPath[i - 1].mIsect.mIsDelta : false;
		if (!lightPath[i].mIsect.mIsDelta && !deltalightVertex)
			sumRi += ri;
	}

	return 1 / (1 + sumRi);
}

real BidirectionalRenderer::PathPdf(
	const Scene						&scene,
	StateSequence					&rand,
	const std::vector<PathVertex>	&path, 
	int								s, 
	int								t)
{
	real p = 1.0;
	for (int i = 0; i < s; ++i) {
		if (i == 0) {
			//real PdfA = Path[0].PdfFwd;
			//const Sphere& light = spheres[numSpheres - 1];
			//real PdfA = 1.0 / (4 * PI * light.rad * light.rad);
			//p *= PdfA;

			real pdfLight, pdfA;
			Light* pLight = scene.SampleOneLight(&pdfLight, rand());
			Intersection isect;
			pLight->SampleOnePoint(&isect, &pdfA, Vec2(rand(), rand()));
			p *= pdfA;
		}
		else if (i == 1) {
			Vec3 lightToHitPoint = (path[1].mIsect.mPos - path[0].mIsect.mPos).Norm();
			real cosThetaLight = std::abs(path[0].mIsect.mNormal.Dot(lightToHitPoint));
			real pdfW = CosineHemispherePdf(cosThetaLight);
			real pdfA = ConvertSolidToArea(pdfW, path[0], path[1]);
			p *= pdfA;
		}
		else {
			if (path[i - 1].mIsect.mIsDelta) p *= 1.0;
			else {
				Vec3 wo = (path[i - 2].mIsect.mPos - path[i - 1].mIsect.mPos).Norm();
				Vec3 wi = (path[i].mIsect.mPos - path[i - 1].mIsect.mPos).Norm();
				real pdfW = path[i - 1].mIsect.mpBSDF->Pdf(wo, wi);
				real pdfA = ConvertSolidToArea(pdfW, path[i - 1], path[i]);
				p *= pdfA;
			}
		}
	}
	if (p == 0.f) return 0;
	for (int i = 0; i < t; ++i) {
		int j = s + t - i - 1;
		if (i == 0) {
			//pinhole
			p *= path[j].mpCamera->PdfPos();
		}
		else if (i == 1) {
			Vec3 cameraToHitPoint = (path[j].mIsect.mPos - path[j + 1].mIsect.mPos).Norm();
			Ray cameraRay(path[j + 1].mIsect.mPos, cameraToHitPoint);
			real pdfW = path[j + 1].mpCamera->PdfDir(cameraRay);
			real pdfA = ConvertSolidToArea(pdfW, path[j + 1], path[j]);
			p *= pdfA;
		}
		else {
			if (path[j + 1].mIsect.mIsDelta) p *= 1.0;
			else {
				Vec3 wo = (path[j + 2].mIsect.mPos - path[j + 1].mIsect.mPos).Norm();
				Vec3 wi = (path[j].mIsect.mPos - path[j + 1].mIsect.mPos).Norm();
				real pdfW = path[j + 1].mIsect.mpBSDF->Pdf(wo, wi);
				real pdfA = ConvertSolidToArea(pdfW, path[j + 1], path[j]);
				p *= pdfA;
			}
		}
	}

	return p;
}

real BidirectionalRenderer::MISWeightV2(
	const Scene					& scene,
	StateSequence				& rand, 
	std::vector<PathVertex>		& lightPath,
	std::vector<PathVertex>		& cameraPath, 
	int							s, 
	int							t,
	PathVertex					& sampled)
{
	if (s + t == 2) return 1.0;

	PathVertex* LightVertex = s > 0 ? &lightPath[s - 1] : nullptr,
		* CameraVertex = t > 0 ? &cameraPath[t - 1] : nullptr;

	ScopedAssignment<PathVertex> a1;
	if (s == 1)
		a1 = { LightVertex, sampled };
	else if (t == 1)
		a1 = { CameraVertex, sampled };

	std::vector<PathVertex> fullPath;
	for (int i = 0; i < s; ++i) fullPath.push_back(lightPath[i]);
	for (int i = t - 1; i >= 0; --i) fullPath.push_back(cameraPath[i]);
	real PdfS = PathPdf(scene, rand, fullPath, s, t);
	real PdfAll = 0;

	for (int nLightVertices = 0; nLightVertices <= s + t - 1; ++nLightVertices) {
		int nCameraVertices = s + t - nLightVertices;
		if (nLightVertices >= 2 && fullPath[nLightVertices - 1].mIsect.mIsDelta) continue;
		if (nLightVertices >= 2 && fullPath[nLightVertices].mIsect.mIsDelta) continue;

		PdfAll += PathPdf(scene, rand, fullPath, nLightVertices, nCameraVertices);

	}
	if ((PdfS == 0.f) || (PdfAll == 0.f)) return 0.f;
	else return std::max(std::min(PdfS / PdfAll, 1.f), 0.f);
}

GYT_NAMESPACE_END

