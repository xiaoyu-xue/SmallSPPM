#pragma once

#include <utility>
#include <cassert>
#include "lowdiscrepency.h"

NAMESPACE_BEGIN

class SamplerEnum {
public:
	//SamplerEnum(){}
	//SamplerEnum(uint32 width, uint32 height){}
	virtual uint64 GetIndex(uint32 sampleNum, uint32 x, uint32 y) const = 0;
	virtual real SampleX(uint32 x, real u) const = 0;
	virtual real SampleY(uint32 y, real v) const = 0;
};


// Determine the index of the i-th sample falling into a pixel, based on the
// elementary interval property of the Halton sequence.
// This is an implementation of the two-dimensional case of the more general
// construction in L. Gruenschloss, M. Raab, and A. Keller: "Enumerating Quasi-Monte
// Carlo Point Sequences in Elementary Intervals".
// This assumes that identity digit permutations are used for the first two components,
// i.e. basis 2 and 3.
class HaltonEnum : public SamplerEnum
{
public:
	// Initialize the enumeration for the given resolution.
	HaltonEnum(uint32 width, uint32 height);

	// Return how many samples per pixel can be queried before sample index overflow occurs.
	uint64 get_max_samples_per_pixel() const { return ~0ull / m_increment; }

	// Return the index of the i-th sample falling into the given pixel (x, y) within the
	// previously given resolution bounds. i must be smaller than the value returned by
	// get_max_samples_per_pixel.
	uint64 get_index(uint64 i, uint32 x, uint32 y) const;

	// Scale the x-component of a sample in [0,1) to [0,width).
	real scale_x(real x) const;

	// Scale the y-component of a sample in [0,1) to [0,height).
	real scale_y(real y) const;

	uint64 GetIndex(uint32 sampleNum, uint32 x, uint32 y) const override {
		return get_index(sampleNum, x, y);
	}

	real SampleX(uint32 x, real u) const override {
		return scale_x(u) - (real)x;
	}

	real SampleY(uint32 y, real v) const override {
		return scale_y(v) - (real)y;
	}

private:
	static std::pair<int, int> extended_euclid(int a, int b);
	static uint64 halton2_inverse(uint64 i, uint32 digits);
	static uint64 halton3_inverse(uint64 i, uint32 digits);

	uint32 m_p2; // Smallest integer with 2^m_p2 >= width.
	uint32 m_p3; // Smallest integer with 3^m_p3 >= height.
	uint32 m_x; // 3^m_p3 * ((2^m_p2)^(-1) mod 3^m_p3).
	uint32 m_y; // 2^m_p2 * ((3^m_p3)^(-1) mod 2^m_p2).
	real m_scale_x; // 2^m_p2.
	real m_scale_y; // 3^m_p3.
	uint64 m_increment; // Product of prime powers, i.e. m_res2 * m_res3.
};

inline HaltonEnum::HaltonEnum(uint32 width, uint32 height)
{
	assert(width && height);

	m_p2 = 0;
	uint32 w = 1;
	while (w < width) // Find 2^m_p2 >= width.
	{
		++m_p2;
		w *= 2;
	}
	m_scale_x = real(w);

	m_p3 = 0;
	uint32 h = 1;
	while (h < height) // Find 3^m_p3 >= height.
	{
		++m_p3;
		h *= 3;
	}
	m_scale_y = real(h);

	m_increment = w * h; // There's exactly one sample per pixel.

	// Determine the multiplicative inverses.
	const std::pair<int, int> inv = extended_euclid(static_cast<int>(h), static_cast<int>(w));
	const uint32 inv2 = (inv.first < 0) ? (inv.first + w) : (inv.first % w);
	const uint32 inv3 = (inv.second < 0) ? (inv.second + h) : (inv.second % h);
	assert((!inv2 && w == 1) || (inv2 > 0 && (h * inv2) % w == 1));
	assert((!inv3 && h == 1) || (inv3 > 0 && (w * inv3) % h == 1));
	m_x = h * inv2;
	m_y = w * inv3;
}

inline uint64 HaltonEnum::get_index(const uint64 i, const uint32 x, const uint32 y) const
{
	// Promote to 64 bits to avoid overflow.
	const uint64 hx = halton2_inverse(x, m_p2);
	const uint64 hy = halton3_inverse(y, m_p3);
	// Apply Chinese remainder theorem.
	const uint64  offset = (hx * m_x + hy * m_y) % m_increment;
	return offset + i * m_increment;
}

inline real HaltonEnum::scale_x(const real x) const
{
	return x * m_scale_x;
}

inline real HaltonEnum::scale_y(const real y) const
{
	return y * m_scale_y;
}

inline std::pair<int, int> HaltonEnum::extended_euclid(const int a, const int b)
{
	if (!b)
		return std::make_pair(1u, 0u);
	const int q = a / b;
	const int r = a % b;
	const std::pair<int, int> st = extended_euclid(b, r);
	return std::make_pair(st.second, st.first - q * st.second);
}

inline uint64 HaltonEnum::halton2_inverse(uint64 index, const uint32 digits)
{
	index = reverse_bit64(index);
	return index >> (64 - digits);
}

inline uint64 HaltonEnum::halton3_inverse(uint64 index, const uint32 digits)
{
	uint32 result = 0;
	for (uint32 d = 0; d < digits; ++d)
	{
		result = result * 3 + index % 3;
		index /= 3;
	}
	return result;
}

NAMESPACE_END