#include "Visibility.h"
#include "Scene.h"
#include "Medium.h"

GY_NAMESPACE_BEGIN

bool VisibilityTester::Unoccluded(const Scene& scene) const {
	Ray ray = p0.SpawnTo(p1);
	return !scene.Intersect(ray);
}


Vec3 VisibilityTester::Tr(const Scene& scene, StateSequence& rand) const {
    Ray ray(p0.SpawnTo(p1));
    Vec3 Tr(1.f);
    while (true) {
        Intersection isect;

        bool hitSurface = scene.Intersect(ray, &isect);
        if (hitSurface && isect.primitive->GetMaterial() != nullptr)
            return Vec3();

        //DEBUG_PIXEL_IF(ThreadIndex()) {
        //    std::cout << "Hit: "<< isect.hit 
        //        << ", Tr: " << Tr << std::endl;
        //}

        if (ray.mpMedium) Tr *= ray.mpMedium->Tr(ray, rand);

        if (!hitSurface) break;
        scene.QueryIntersectionInfo(ray, &isect);
        ray = isect.SpawnTo(p1);
    }
    return Tr;
}

GY_NAMESPACE_END