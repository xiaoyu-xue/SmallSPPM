#pragma once

#include "BidirectionalRenderer.h"
#include "core/Scene.h"
#include "core/Intersection.h"
#include "light/Light.h"

GYT_NAMESPACE_BEGIN

int Gyunity::BidirectionalRenderer::GenerateLightPath(const Scene& scene, SimpleSampler& sampler, std::vector<PathVertex>& lightPath, int maxDepth)
{
	real pdfLight;
	Light* pLight = scene.SampleOneLight(&pdfLight, sampler.Get1D());
	Intersection isect;
	Vec3 dir;
	real pdfPos, pdfDir;
	Vec3 Le = pLight->Emission();
	pLight->SampleOnLight(&isect, &dir, &pdfPos, &pdfDir, sampler.Get2D(), sampler.Get2D());
	Ray ray(isect.mPos + dir * RayEps, dir);
	lightPath[0].mIsect.mPos = isect.mPos;
	lightPath[0].mIsect.mNormal = isect.mNormal;
	lightPath[0].mThroughput = Le;
	lightPath[0].mPdfFwd = pdfPos;
	lightPath[0].mIsect.mIsDelta = false;
	lightPath[0].mpLight = pLight;
	real cosTheta = dir.Norm().Dot(isect.mNormal);
	Vec3 throughPut = lightPath[0].mThroughput * std::abs(cosTheta) / pdfPos / pdfDir / pdfLight;
}


real BidirectionalRenderer::ConvertSolidToArea(real pdfW, const PathVertex& Vertex, const PathVertex& nxt)
{
	Vec3 dir = nxt.mIsect.mPos - Vertex.mIsect.mPos;
	real dist2 = dir.Length2();
	if (dist2 == 0) return 0;
	dir.Normalize();
	real cosTheta = std::abs(dir.Dot(nxt.mIsect.mNormal));
	return pdfW * cosTheta / dist2;
}


GYT_NAMESPACE_END

