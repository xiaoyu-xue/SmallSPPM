#pragma once

#include <utility>
#include <cassert>
#include "lowdiscrepency.h"

class SamplerEnum {
public:
	//SamplerEnum(){}
	//SamplerEnum(unsigned width, unsigned height){}
	virtual unsigned long long GetIndex(unsigned sampleNum, unsigned x, unsigned y) const = 0;
	virtual double SampleX(unsigned x, double u) const = 0;
	virtual double SampleY(unsigned y, double v) const = 0;
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
	HaltonEnum(unsigned width, unsigned height);

	// Return how many samples per pixel can be queried before sample index overflow occurs.
	unsigned long long get_max_samples_per_pixel() const { return ~0ull / m_increment; }

	// Return the index of the i-th sample falling into the given pixel (x, y) within the
	// previously given resolution bounds. i must be smaller than the value returned by
	// get_max_samples_per_pixel.
	unsigned long long get_index(unsigned long long i, unsigned x, unsigned y) const;

	// Scale the x-component of a sample in [0,1) to [0,width).
	float scale_x(float x) const;

	// Scale the y-component of a sample in [0,1) to [0,height).
	float scale_y(float y) const;

	unsigned long long GetIndex(unsigned sampleNum, unsigned x, unsigned y) const override {
		return get_index(sampleNum, x, y);
	}

	double SampleX(unsigned x, double u) const override {
		return scale_x(u) - (float)x;
	}

	double SampleY(unsigned y, double v) const override {
		return scale_y(v) - (float)y;
	}

private:
	static std::pair<int, int> extended_euclid(int a, int b);
	static unsigned long long halton2_inverse(unsigned long long i, unsigned digits);
	static unsigned long long halton3_inverse(unsigned long long i, unsigned digits);

	unsigned m_p2; // Smallest integer with 2^m_p2 >= width.
	unsigned m_p3; // Smallest integer with 3^m_p3 >= height.
	unsigned m_x; // 3^m_p3 * ((2^m_p2)^(-1) mod 3^m_p3).
	unsigned m_y; // 2^m_p2 * ((3^m_p3)^(-1) mod 2^m_p2).
	float m_scale_x; // 2^m_p2.
	float m_scale_y; // 3^m_p3.
	unsigned long long m_increment; // Product of prime powers, i.e. m_res2 * m_res3.
};

HaltonEnum::HaltonEnum(unsigned width, unsigned height)
{
	assert(width && height);

	m_p2 = 0;
	unsigned w = 1;
	while (w < width) // Find 2^m_p2 >= width.
	{
		++m_p2;
		w *= 2;
	}
	m_scale_x = float(w);

	m_p3 = 0;
	unsigned h = 1;
	while (h < height) // Find 3^m_p3 >= height.
	{
		++m_p3;
		h *= 3;
	}
	m_scale_y = float(h);

	m_increment = w * h; // There's exactly one sample per pixel.

	// Determine the multiplicative inverses.
	const std::pair<int, int> inv = extended_euclid(static_cast<int>(h), static_cast<int>(w));
	const unsigned inv2 = (inv.first < 0) ? (inv.first + w) : (inv.first % w);
	const unsigned inv3 = (inv.second < 0) ? (inv.second + h) : (inv.second % h);
	assert((!inv2 && w == 1) || (inv2 > 0 && (h * inv2) % w == 1));
	assert((!inv3 && h == 1) || (inv3 > 0 && (w * inv3) % h == 1));
	m_x = h * inv2;
	m_y = w * inv3;
}

inline unsigned long long HaltonEnum::get_index(const unsigned long long i, const unsigned x, const unsigned y) const
{
	// Promote to 64 bits to avoid overflow.
	const unsigned long long hx = halton2_inverse(x, m_p2);
	const unsigned long long hy = halton3_inverse(y, m_p3);
	// Apply Chinese remainder theorem.
	const unsigned long long  offset = (hx * m_x + hy * m_y) % m_increment;
	return offset + i * m_increment;
}

inline float HaltonEnum::scale_x(const float x) const
{
	return x * m_scale_x;
}

inline float HaltonEnum::scale_y(const float y) const
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

inline unsigned long long HaltonEnum::halton2_inverse(unsigned long long index, const unsigned digits)
{
	index = reverse_bit64(index);
	return index >> (64 - digits);
}

inline unsigned long long HaltonEnum::halton3_inverse(unsigned long long index, const unsigned digits)
{
	unsigned result = 0;
	for (unsigned d = 0; d < digits; ++d)
	{
		result = result * 3 + index % 3;
		index /= 3;
	}
	return result;
}
