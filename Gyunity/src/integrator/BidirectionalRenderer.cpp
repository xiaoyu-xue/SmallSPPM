#pragma once

#include "BidirectionalRenderer.h"
#include "core/Scene.h"
#include "core/Intersection.h"
#include "light/Light.h"

GYT_NAMESPACE_BEGIN

int Gyunity::BidirectionalRenderer::GenerateLightPath(const Scene& scene, StateSequence& rand, MemoryPool& arena, std::vector<PathVertex>&lightPath, int maxDepth)
{
	real pdfLight;
	Light* pLight = scene.SampleOneLight(&pdfLight, rand());
	Intersection isect;
	Vec3 dir;
	real pdfPos, pdfDir;
	Vec3 Le = pLight->Emission();
	pLight->SampleOnLight(&isect, &dir, &pdfPos, &pdfDir, Vec2(rand(), rand()), Vec2(rand(), rand()));
	lightPath[0].mIsect.mPos = isect.mPos;
	lightPath[0].mIsect.mNormal = isect.mNormal;
	lightPath[0].mThroughput = Le;
	lightPath[0].mPdfFwd = pdfPos;
	lightPath[0].mIsect.mIsDelta = false;
	lightPath[0].mpLight = pLight;
	real cosTheta = dir.Norm().Dot(isect.mNormal);
	Vec3 throughPut = Le * std::abs(cosTheta) / pdfPos / pdfDir / pdfLight;
	Ray ray(isect.mPos + dir * RayEps, dir);
	//int nLightVertices = BidirectionalRenderer::Trace(scene, arena, rand, ray, throughPut, pdfDir, lightPath, 1, maxDepth, TransportMode::Importance);
	int nLightVertices = BidirectionalRenderer::TraceV2(scene, arena, rand, ray, 1, throughPut, pdfDir, lightPath, maxDepth, TransportMode::Importance);
	return nLightVertices;
}


int BidirectionalRenderer::GenerateCameraPath(const Camera& camera, StateSequence& rand, MemoryPool& arena, std::vector<PathVertex>& cameraPath, const Ray& cameraRay, int maxdDpth)
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
	
	return 1;
}

int BidirectionalRenderer::Trace(
	const Scene				&scene,
	MemoryPool				&arena,
	StateSequence			&rand,
	const Ray				&ray,
	Vec3					throughput,
	real					pdfFwd,
	std::vector<PathVertex>	&path,
	int						depth,
	int						maxDepth,
	TransportMode			mode)
{
	Ray r = ray;
	int bounce = depth;
	real pdfW = pdfFwd;
	while (1) {
		PathVertex& prev = path[bounce - 1];
		PathVertex& vertex = path[bounce];
		Intersection& isect = path[bounce].mIsect;

		if (!scene.Intersect(r, &isect))
		{
			return bounce - 1;
		}
		scene.QueryIntersectionInfo(r, &isect);
		isect.ComputeScatteringFunction(arena, mode);

		path[bounce].mIsDelta = isect.mpBSDF->IsDelta();
		path[bounce].mThroughput = throughput;
		path[bounce].mPdfFwd = ConvertSolidToArea(pdfW, prev, vertex);
		//path[bounce].mPdfFwd = pdfW;
		++bounce;
		if (bounce >= maxDepth) break;


		Vec3 wo;
		Vec3 f = isect.mpBSDF->Sample(-1 * r.mDir, &wo, &pdfW, Vec3(rand(), rand(), rand()));
		wo.Normalize();
		throughput = throughput * f * (std::abs(wo.Dot(isect.mNormal))) / pdfW;

		real pdfWPrev = isect.mpBSDF->Pdf(wo, -1 * r.mDir);
		prev.mPdfPrev = ConvertSolidToArea(pdfWPrev, vertex, prev);

		if (isect.mIsDelta) {
			pdfW = pdfWPrev = 0;
		}

		r = Ray(isect.mPos + wo * RayEps, wo);
	}
	return bounce - 1;
}

int BidirectionalRenderer::TraceV2(
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
		path[bounce].mPdfFwd = ConvertSolidToArea(pdfDir, prev, vertex);;
		path[bounce].mThroughput = throughput;

		bounce++;
		if (bounce >= maxDepth) break;

		Vec3 wo;
		Vec3 f = isect.mpBSDF->Sample(-1 * ray.mDir, &wo, &pdfDir, Vec3(rand(), rand(), rand()));
		throughput = throughput * f * std::abs(Dot(wo, isect.mNormal)) / pdfDir;

		real pdfWPrev = isect.mpBSDF->Pdf(wo, -1 * r.mDir);
		prev.mPdfPrev = ConvertSolidToArea(pdfWPrev, vertex, prev);

		if (isect.mIsDelta) {
			pdfDir = pdfWPrev = 0;
		}

		ray = Ray(isect.mPos + wo * RayEps, wo);


	}
	return bounce - 1;
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

