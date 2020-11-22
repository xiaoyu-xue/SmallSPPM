#pragma once
#include "system/memory.h"

NAMESPACE_BEGIN

enum ReflectionType { DIFF, SPEC, REFR };  // material types, used in radiance()
enum class TransportMode { Radiance = 1, Importance };

class Intersection;

class Material {
public:
	virtual ~Material() { }
	//virtual void ComputeScatteringFunction(Intersection* isect,
	//	TransportMode mode = TransportMode::Radiance) const = 0;
	virtual void ComputeScatteringFunction(Intersection *isect, MemoryArena& arena,
		TransportMode mode = TransportMode::Radiance) const = 0;
};

NAMESPACE_END
