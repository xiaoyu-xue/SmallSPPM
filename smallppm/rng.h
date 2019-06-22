#pragma once
#include <random>
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
private:
	std::mt19937_64 rng;
	std::uniform_real_distribution<double> uniform;
};