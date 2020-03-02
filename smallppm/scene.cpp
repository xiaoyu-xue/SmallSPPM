#include "scene.h"

NAMESPACE_BEGIN

void Scene::QueryIntersectionInfo(const Ray& ray, Intersection* isect) const {
	primitives[isect->shapeId]->QueryIntersectionInfo(ray, isect);
}

void Scene::GetBoundingSphere(Vec3* center, real* radius) const {
	accelerator->WorldBound().GetBoundingSphere(center, radius);
}


NAMESPACE_BEGIN