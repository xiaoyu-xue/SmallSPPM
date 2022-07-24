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
	static int GenerateLightPath(
		const Scene				&scene, 
		StateSequence			&rand, 
		MemoryPool				&arena, 
		std::vector<PathVertex>	&lightPath, 
		int						maxdepth);

	static int GenerateCameraPath(
		const Scene				&scene,
		const Camera			&camera, 
		StateSequence			&rand, 
		MemoryPool				&arena, 
		std::vector<PathVertex>	&cameraPath, 
		const Ray				&cameraRay, 
		int						maxdepth);

	static int Trace(
		const Scene				&scene,
		MemoryPool				&arena,
		StateSequence			&rand,
		const Ray				&r,
		int						depth,
		Vec3					throughput,
		real					pdfFwd,
		std::vector<PathVertex>	&lightPath,
		int						maxDepth,
		TransportMode			mode);

	static real ConvertSolidToArea(
		real				pdfW, 
		const PathVertex	&vertex, 
		const PathVertex	&nxt);

	static bool IsConnectable(
		const Scene		&scene, 
		const Vec3		&pointA, 
		const Vec3		&pointB);

	static real G(
		const PathVertex& vertexA, 
		const PathVertex& vertexB);

	static real MISWeight(
		const Scene				&scene,
		StateSequence			&rand,
		std::vector<PathVertex>	&lightPath, 
		std::vector<PathVertex>	&caameraPath,
		int						s, 
		int						t, 
		PathVertex				&sampled);

	static real MISWeightV2(
		const Scene				&scene, 
		StateSequence			&rand, 
		std::vector<PathVertex>	&lightPath, 
		std::vector<PathVertex>	&cameraPath, 
		int						s, 
		int						t, 
		PathVertex				&sampled);

	static real PathPdf(
		const Scene						&scene, 
		StateSequence					&rand, 
		const std::vector<PathVertex>	&path, 
		int								s, 
		int								t);
};

GYT_NAMESPACE_END