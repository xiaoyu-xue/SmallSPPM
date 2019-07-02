#pragma once

#include <memory>
#include <random>
#include "lowdiscrepency.h"

NAMESPACE_BEGIN

class Sampler {
public:
	virtual real Sample(unsigned int d, unsigned long long i) = 0;
};

class StateSequence {
protected:
	int cursor = 0;

public:
	virtual real Sample() = 0;

	virtual real operator()() {
		return Sample();
	}

	int GetCursor() const {
		return cursor;
	}
	/*
	void assert_cursor_pos(int cursor) const {
		assert_info(
			this->cursor == cursor,
			std::string("Cursor position should be " + std::to_string(cursor) +
				" instead of " + std::to_string(this->cursor)));
	}*/



};

class RandomStateSequence : public StateSequence {
private:
	std::shared_ptr<Sampler> sampler;
	long long instance;

public:
	RandomStateSequence(std::shared_ptr<Sampler> sampler, long long instance)
		: sampler(sampler), instance(instance) {
	}

	real Sample() override {
		//assert_info(sampler != nullptr, "null sampler");
		real ret = sampler->Sample(cursor++, instance);
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
	RandomSampler(int seed) : rng(seed) {

	}

	real Sample(unsigned int d, unsigned long long i) {
		return uniform(rng);
		//return Rand();
	}
private:
	std::mt19937_64 rng;
	std::uniform_real_distribution<real> uniform;
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


NAMESPACE_END