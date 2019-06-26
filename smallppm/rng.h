#pragma once
#include <random>

NAMESPACE_BEGIN

class Rng {
public:
	Rng(int seed = 123) : rng(seed){}

	double operator()() {
		return uniform(rng);
	}

	double operator()(unsigned long long i) {
		return uniform(rng);
	}

	double GetDouble() {
		return uniform(rng);
	}

	std::mt19937_64 GetRngEngine() {
		return rng;
	}

private:
	std::mt19937_64 rng;
	std::uniform_real_distribution<double> uniform;
};

NAMESPACE_END