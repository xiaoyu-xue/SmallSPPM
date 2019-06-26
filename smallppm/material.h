#pragma once

NAMESPACE_BEGIN

enum ReflectionType { DIFF, SPEC, REFR };  // material types, used in radiance()
enum class TransportMode { Radiance = 1, Importance };

NAMESPACE_END
