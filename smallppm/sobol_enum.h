#pragma once
#include "sampler_enum.h"
#include "scalar.h"

NAMESPACE_BEGIN

extern const unsigned long long vdc_sobol_matrices[][52];
extern const unsigned long long vdc_sobol_matrices_inv[][52];

class SobolEnum : public SamplerEnum {
public:
	SobolEnum() {}

	SobolEnum(uint32 width, uint32 height) : resX(width), resY(height) {
		mResolution = RoundUpPow2(std::max((int32)resX,  (int32)resY));
		mLog2Resolution = Log2Int(mResolution);
		assert(1 << mLog2Resolution == mResolution);
	}

	uint64 GetIndex(uint32 sampleNum, uint32 x, uint32 y) const {
        const uint64 scramble_x = 0;
        const uint64 scramble_y = 0;
        uint64 index = look_up(mLog2Resolution, sampleNum, x, y, scramble_x, scramble_y);
        return index;
	}

	real SampleX(uint32 x, real u) const {
		return mResolution * u - real(x);
	}

	real SampleY(uint32 y, real v) const {
		return mResolution * v - real(y);
	}

private:
    // Return the index of the frame-th sample falling
    // into the square elementary interval (px, py),
    // without using look-up tables.
    inline uint64 look_up(
        const uint32 m,
        uint32 frame,
        uint32 px,
        uint32 py,
        const uint64 scramble_x,
        const uint64 scramble_y) const
    {
        if (m == 0) return 0;

        const uint32 m2 = m << 1;
        uint64 index = uint64(frame) << m2;

        // Note: the delta value only depends on frame
        // and m, thus it can be cached across multiple
        // function calls, if desired.
        uint64 delta = 0;
        for (uint32 c = 0; frame; frame >>= 1, ++c)
            if (frame & 1) // Add flipped column m + c + 1.
                delta ^= vdc_sobol_matrices[m - 1][c];

        px ^= scramble_x >> (64 - m);
        py ^= scramble_y >> (64 - m);
        uint64 b = ((uint64(px) << m) | py) ^ delta; // flipped b

        for (uint32 c = 0; b; b >>= 1, ++c)
            if (b & 1) // Add column 2 * m - c.
                index ^= vdc_sobol_matrices_inv[m - 1][c];

        return index;
    }

	uint32 resX, resY;
	int mResolution, mLog2Resolution;
};

NAMESPACE_END