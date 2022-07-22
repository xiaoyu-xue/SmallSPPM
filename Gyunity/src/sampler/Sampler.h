#pragma once

#include <memory>
#include <random>
#include "math/Linagl.h"
#include "math/LowDiscrepency.h"
#include "Rng.h"

GYT_NAMESPACE_BEGIN

class Sampler {
public:
	virtual real Sample(uint32 d, uint64 i) = 0;
	virtual std::shared_ptr<Sampler> Clone(uint32 seed) = 0;
};

class SimpleSampler {
public:
	SimpleSampler(int seed = 1234) : rng(seed) {}

	real Get1D() {
		return rng.GetFloat();
	}

	Vec2 Get2D() {
		return Vec2(Get1D(), Get1D());
	}

	Vec3 Get3D() {
		return Vec3(Get1D(), Get1D(), Get1D());
	}

	std::shared_ptr<SimpleSampler> Clone(uint32 seed) {
		return std::make_shared<SimpleSampler>(seed);
	}

private:
	Rng rng;
};

class StateSequence {
protected:
	int mCursor = 0;

public:
	virtual real Sample() = 0;

	virtual real operator()() {
		return Sample();
	}

	int GetCursor() const {
		return mCursor;
	}
};

class RandomStateSequence : public StateSequence {
private:
	std::shared_ptr<Sampler> mpSampler;
	int64 mInstance;

public:
	RandomStateSequence(std::shared_ptr<Sampler> sampler, int64 instance)
		: mpSampler(sampler), mInstance(instance) 
	{
	}

	real Sample() override {
		//assert_info(sampler != nullptr, "null sampler");
		real ret = mpSampler->Sample(mCursor++, mInstance);
		//assert_info(ret >= 0, "sampler output should be non-neg");
		if (ret > 1 + 1e-5f) {
			printf("Warning: sampler returns value > 1: [%f]", ret);
		}
		if (ret >= 1) {
			ret = 0;
		}
		return ret;
	}

};

class RandomSampler : public Sampler {
public:
	RandomSampler(int seed) : rng(seed), pbrtRng(seed) {

	}

	real Sample(unsigned int d, unsigned long long i) {
		return uniform(rng);
		//return pbrtRng.UniformFloat();
	}

	std::shared_ptr<Sampler> Clone(uint32 seed) override {
		RandomSampler* newSampler = new RandomSampler(seed);
		newSampler->pbrtRng.SetSequence(seed);
		return std::shared_ptr<Sampler>(newSampler);
	}
private:
	std::mt19937_64 rng;
	std::uniform_real_distribution<real> uniform;
	RNG pbrtRng;
};

class PrimeList {
public:
	PrimeList() {
		for (int i = 2; i <= 10000; i++) {
			if (IsPrime(i)) {
				primes.push_back(i);
			}
		}
		ASSERT(primes.size() == 1229);
	}

	int GetPrime(int i) {
		return primes[i];
	}

	int GetPrimesNum() {
		return (int)primes.size();
	}

private:
	std::vector<int> primes;
};


class RegularHaltonSampler : public Sampler {
public:
	real Sample(unsigned int d, unsigned long long i) override {
		ASSERT(d < (unsigned)primeList.GetPrimesNum());
		real val = hal(d, i + 1);  // The first one is evil...
		return val;
	}

	std::shared_ptr<Sampler> Clone(uint32 seed) override {
		return std::shared_ptr<Sampler>(new RegularHaltonSampler(*this));
	}

private:
	inline int rev(const int i, const int p) const {
		return i == 0 ? i : p - i;
	}

	real hal(const int d, long long j) const {
		const int p = primeList.GetPrime(d);
		real h = 0.f, f = 1.f / p, fct = f;
		while (j > 0) {
			h += rev(j % p, p) * fct;
			j /= p;
			fct *= f;
		}
		return h;
	}
	static PrimeList primeList;
};


GYT_NAMESPACE_END