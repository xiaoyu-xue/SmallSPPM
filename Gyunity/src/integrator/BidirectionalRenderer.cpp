#pragma once

#include "BidirectionalRenderer.h"
#include "core/Scene.h"
#include "core/Intersection.h"
#include "light/Light.h"

GYT_NAMESPACE_BEGIN

int Gyunity::BidirectionalRenderer::GenerateLightPath(
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


int BidirectionalRenderer::GenerateCameraPath(const Scene &scene, const Camera& camera, StateSequence& rand, MemoryPool& arena, std::vector<PathVertex>& cameraPath, const Ray& cameraRay, int maxdDpth)
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
	Ray ray = r;
	int bounce = depth;
	real pdfDir = pdfFwd;
	while (true) {
		PathVertex& prev = path[bounce - 1];
		PathVertex& vertex = path[bounce];
		Intersection& isect = path[bounce].mIsect;

		if (!scene.Intersect(ray, &isect)) break;

		scene.QueryIntersectionInfo(ray, &isect);
		isect.ComputeScatteringFunction(arena, mode);

		path[bounce].mIsDelta = isect.mpBSDF->IsDelta();
		path[bounce].mPdfFwd = BidirectionalRenderer::ConvertSolidToArea(pdfDir, prev, vertex);
		path[bounce].mThroughput = throughput;

		bounce++;
		if (bounce >= maxDepth) break;

		Vec3 wi;
		Vec3 f = isect.mpBSDF->Sample(-1 * ray.mDir, &wi, &pdfDir, Vec3(rand(), rand(), rand()));
		
		throughput = throughput * f * std::abs(Dot(wi, isect.mNormal)) / pdfDir;
		real pdfWPrev = isect.mpBSDF->Pdf(wi, -1 * ray.mDir);
		prev.mPdfPrev = BidirectionalRenderer::ConvertSolidToArea(pdfWPrev, vertex, prev);

		ray = Ray(isect.mPos + wi * RayEps, wi);
	}
	return bounce;
}

real BidirectionalRenderer::ConvertSolidToArea(real pdfW, const PathVertex& vertex, const PathVertex& nxt)
{
	Vec3 dir = nxt.mIsect.mPos - vertex.mIsect.mPos;
	real dist2 = dir.Length2();
	if (dist2 == 0) return 0;
	dir.Normalize();
	real cosTheta = std::abs(dir.Dot(nxt.mIsect.mNormal));
	return pdfW * cosTheta / dist2;
}

bool BidirectionalRenderer::IsConnectable(const Scene &scene, const Vec3& pointA, const Vec3& pointB)
{
	Vec3 dir = pointA - pointB;
	real dist = dir.Length();
	Ray ray(pointB, dir.Norm(), 0.0, dist - Eps);  //Be careful ! Float error will cause the incorrect intersection

#ifdef _DEBUG
	if (scene.Intersect(ray)) return false;
#else
	if (scene.Intersect(ray)) return false;
#endif
	return true;
}



GYT_NAMESPACE_END

