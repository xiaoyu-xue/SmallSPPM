#pragma once

#include <utility>
#include <cassert>
#include "lowdiscrepency.h"

NAMESPACE_BEGIN


extern const unsigned long long vdc_sobol_matrices[][52];
extern const unsigned long long vdc_sobol_matrices_inv[][52];


class SamplerEnum {
public:
	SamplerEnum(){}
	SamplerEnum(uint32 width, uint32 height): resX(width), resY(height) {}
	virtual uint64 GetIndex(uint32 sampleNum, uint32 x, uint32 y) const {
		return sampleNum * y;
	}
	virtual real SampleX(uint32 x, real u) const {
		return u;
	}
	virtual real SampleY(uint32 y, real v) const {
		return v;
	}
private:
	uint32 resX, resY;
};

NAMESPACE_END