#pragma once
#include "system/Memory.h"

GYT_NAMESPACE_BEGIN

enum ReflectionType { DIFF, SPEC, REFR };  // material types, used in radiance()
enum class TransportMode { Radiance = 1, Importance = 2 };

class Intersection;

class Material {
public:
	virtual ~Material() { }
	//virtual void ComputeScatteringFunction(Intersection* isect,
	//	TransportMode mode = TransportMode::Radiance) const = 0;
	virtual void ComputeScatteringFunction(Intersection *isect, MemoryPool& arena,
		TransportMode mode = TransportMode::Radiance) const = 0;
};

GYT_NAMESPACE_END
