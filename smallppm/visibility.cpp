#include "visibility.h"
#include "scene.h"
#include "medium.h"

NAMESPACE_BEGIN

bool VisibilityTester::Unoccluded(const Scene& scene) const {
	Ray ray = p0.SpawnTo(p1);
	Intersection isect;
	return scene.Intersect(ray);
}


Vec3 VisibilityTester::Tr(const Scene& scene, StateSequence& rand) const {
    Ray ray(p0.SpawnTo(p1));
    Vec3 Tr(1.f);
    while (true) {
        Intersection isect;

        bool hitSurface = scene.Intersect(ray, &isect);
        if (hitSurface && isect.primitive->GetMaterial() != nullptr)
            return Vec3();

        if (ray.medium) Tr *= ray.medium->Tr(ray, rand);

        if (!hitSurface) break;
        ray = isect.SpawnTo(p1);
    }
    return Tr;
}

NAMESPACE_END