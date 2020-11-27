#include "Scene.h"

NAMESPACE_BEGIN

void Scene::QueryIntersectionInfo(const Ray& ray, Intersection* isect) const {
	primitives[isect->primId]->QueryIntersectionInfo(ray, isect);
}

void Scene::GetBoundingSphere(Vec3* center, real* radius) const {
	accelerator->WorldBound().GetBoundingSphere(center, radius);
}


bool Scene::IntersectTr(Ray& ray, StateSequence& rand, Intersection* isect, Vec3* Tr) const {
    *Tr = Vec3(1.f);
    while (true) {
        bool hitSurface = Intersect(ray, isect);
        // Accumulate beam transmittance for ray segment
        if (ray.mpMedium) *Tr *= ray.mpMedium->Tr(ray, rand);

        // Initialize next ray segment or terminate transmittance computation
        if (!hitSurface) return false;
        if (isect->primitive->GetMaterial() != nullptr) return true;
        ray = isect->SpawnRay(ray.mDir);
    }
}

NAMESPACE_END