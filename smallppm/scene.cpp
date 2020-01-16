#include "scene.h"

NAMESPACE_BEGIN

void Scene::QueryIntersectionInfo(const Ray& ray, Intersection* isect) const {
	primitives[isect->shapeId]->QueryIntersectionInfo(ray, isect);
}

NAMESPACE_BEGIN