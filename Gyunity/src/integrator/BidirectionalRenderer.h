#pragma once

#include "ForwardDecl.h"
#include "PathVertex.h"

GYT_NAMESPACE_BEGIN

template <typename Type>
class ScopedAssignment {
public:
	// ScopedAssignment Public Methods
	ScopedAssignment(Type* target = nullptr, Type value = Type())
		: target(target) {
		if (target) {
			backup = *target;
			*target = value;
		}
	}
	~ScopedAssignment() {
		if (target) *target = backup;
	}
	ScopedAssignment(const ScopedAssignment&) = delete;
	ScopedAssignment& operator=(const ScopedAssignment&) = delete;
	ScopedAssignment& operator=(ScopedAssignment&& other) {
		target = other.target;
		backup = other.backup;
		other.target = nullptr;
		return *this;
	}

private:
	Type* target, backup;
};

class BidirectionalRenderer {
public:
	static int GenerateLightPath(const Scene& scene, SimpleSampler& sampler, std::vector<PathVertex>& LightPath, int maxdepth);

	static int GenerateCameraPath(const Camera& camera, Sampler& sampler, std::vector<PathVertex>& CameraPath, const Ray& cameraRay, int maxdepth);

	real ConvertSolidToArea(real pdfW, const PathVertex& Vertex, const PathVertex& nxt);
};

GYT_NAMESPACE_END