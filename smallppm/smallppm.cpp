#include "svpng.inc"
#include <math.h>  
#include <stdlib.h> 
#include <stdio.h>  
#include <random>
#include <vector>
#include <crtdbg.h>
#include <iostream>
#include <algorithm>
#include <atomic>
#include <mutex>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string.h>
#include "halton_sampler.h"
#define _CRTDBG_MAP_ALLOC

#define ASSERT(expr) \
	do { if(!(expr)) { std::cerr << "Error: assertion `"#expr"' failed at " << __FILE__ << ":" << __LINE__ << std::endl; exit(2); } } while(0)

const double PI = 3.14159265358979;
const double INV_PI = 0.31830988618379067154;
const double ALPHA = 0.66666667;
const long long  render_stage_number = 700000;
const double PiOver2 = 1.57079632679489661923;
const double PiOver4 = 0.78539816339744830961;
const double eps = 1e-6;
const double Inf = 1e20;
const double rayeps = 1e-3;

class Shape;
class BSDF;

/*
std::mt19937_64 rng(1234);
std::uniform_real_distribution<double> uniform;
double Random() {
	return uniform(rng);
}
*/

enum ReflectionType { DIFF, SPEC, REFR };  // material types, used in radiance()

int IsPrime(int a) noexcept {
	ASSERT(a >= 2);
	for (int i = 2; i * i <= a; i++) {
		if (a % i == 0)
			return false;
	}
	return true;
}

// Halton sequence with reverse permutation
/*
int primes[61] = {
	2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,
	83,89,97,101,103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,
	191,193,197,199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283
};
inline int rev(const int i, const int p) {
	if (i == 0) return i; else return p - i;
}

double hal(const int b, int j) {
	const int p = primes[b];
	double h = 0.0, f = 1.0 / (double)p, fct = f;
	while (j > 0) {
		h += rev(j % p, p) * fct; j /= p; fct *= f;
	}
	return h;
}*/

template <class T1, class T2, class T3>
T1 Clamp(const T1& tVal, const T2& tMin, const T3& max)
{
	if (tVal < tMin) return tMin;
	if (tVal > max) return max;
	return tVal;
}

class Spinlock {
protected:
	std::atomic<bool> latch;

public:
	Spinlock() : Spinlock(false) {
	}

	Spinlock(bool flag) {
		latch.store(flag);
	}

	Spinlock(int flag) : Spinlock(flag != 0) {
	}

	void lock() {
		bool unlatched = false;
		while (!latch.compare_exchange_weak(unlatched, true,
			std::memory_order_acquire)) {
			unlatched = false;
		}
	}

	void unlock() {
		latch.store(false, std::memory_order_release);
	}

	Spinlock(const Spinlock &o) {
		// We just ignore racing condition here...
		latch.store(o.latch.load());
	}

	Spinlock &operator=(const Spinlock &o) {
		// We just ignore racing condition here...
		latch.store(o.latch.load());
		return *this;
	}
};

struct Distribution1D
{
public:
	Distribution1D(const double *f, int n)
	{
		count = n;
		func = new double[n];
		memcpy(func, f, n * sizeof(double));
		cdf = new double[n + 1];

		// Compute integral of step function at $x_i$
		cdf[0] = 0.;
		for (int i = 1; i < count + 1; ++i)
			cdf[i] = cdf[i - 1] + func[i - 1] / n;

		// Transform step function integral into CDF
		funcInt = cdf[count];
		if (funcInt == 0.f)
		{
			for (int i = 1; i < n + 1; ++i)
				cdf[i] = double(i) / double(n);
		}
		else
		{
			for (int i = 1; i < n + 1; ++i)
				cdf[i] /= funcInt;
		}
	}

	~Distribution1D()
	{
		delete[] func;
		delete[] cdf;
	}

	double SampleContinuous(double u, double *pdf, int *off = NULL) const
	{
		// Find surrounding CDF segments and _offset_
		double *ptr = std::upper_bound(cdf, cdf + count + 1, u);
		int offset = Clamp(int(ptr - cdf - 1), 0, count - 1);
		if (off) *off = offset;
		ASSERT(offset < count);
		ASSERT(u >= cdf[offset] && (u < cdf[offset + 1] || u == 1));

		// Fix the case when func ends with zeros
		if (cdf[offset] == cdf[offset + 1])
		{
			ASSERT(u == 1.0f);

			do { offset--; } while (cdf[offset] == cdf[offset + 1] && offset > 0);

			ASSERT(cdf[offset] != cdf[offset + 1]);
		}

		// Compute offset along CDF segment
		double du = (u - cdf[offset]) / (cdf[offset + 1] - cdf[offset]);
		ASSERT(!std::isnan(du));
		
		// Compute PDF for sampled offset
		if (pdf) *pdf = func[offset] / funcInt;
		ASSERT(func[offset] > 0);

		// Return $x\in{}[0,1]$ corresponding to sample
		return (offset + du) / count;
	}

	int SampleDiscrete(double u, double *pdf) const
	{
		// Find surrounding CDF segments and _offset_
		double *ptr = std::upper_bound(cdf, cdf + count + 1, u);
		int offset = std::max(0, int(ptr - cdf - 1));
		ASSERT(offset < count);
		ASSERT(u >= cdf[offset] && u < cdf[offset + 1]);
		if (pdf) *pdf = func[offset] / (funcInt * count);
		return offset;
	}

private:
	friend struct Distribution2D;
	double *func, *cdf;
	double funcInt;
	int count;
};

struct Distribution2D
{
public:
	Distribution2D(const double *func, int nu, int nv)
	{
		pConditionalV.reserve(nv);
		for (int v = 0; v < nv; ++v)
		{
			// Compute conditional sampling distribution for $\tilde{v}$
			pConditionalV.push_back(new Distribution1D(&func[v*nu], nu));
		}

		// Compute marginal sampling distribution $p[\tilde{v}]$
		std::vector<double> marginalFunc;
		marginalFunc.reserve(nv);
		for (int v = 0; v < nv; ++v)
			marginalFunc.push_back(pConditionalV[v]->funcInt);
		pMarginal = new Distribution1D(&marginalFunc[0], nv);
	}

	~Distribution2D()
	{
		delete pMarginal;
		for (uint32_t i = 0; i < pConditionalV.size(); ++i)
			delete pConditionalV[i];
	}

	void SampleContinuous(double u0, double u1, double uv[2], double *pdf) const
	{
		double pdfs[2];
		int v;
		uv[1] = pMarginal->SampleContinuous(u1, &pdfs[1], &v);
		uv[0] = pConditionalV[v]->SampleContinuous(u0, &pdfs[0]);
		*pdf = pdfs[0] * pdfs[1];
	}

	double Pdf(double u, double v) const
	{
		int iu = Clamp((int)(u * pConditionalV[0]->count), 0,
			pConditionalV[0]->count - 1);
		int iv = Clamp((int)(v * pMarginal->count), 0,
			pMarginal->count - 1);
		if (pConditionalV[iv]->funcInt * pMarginal->funcInt == 0.f) return 0.f;
		return (pConditionalV[iv]->func[iu] * pMarginal->func[iv]) /
			(pConditionalV[iv]->funcInt * pMarginal->funcInt);
	}

private:
	std::vector<Distribution1D *> pConditionalV;
	Distribution1D *pMarginal;
};


class Sampler {
public:
	virtual double Sample(int d, long long i) = 0;
};

class StateSequence {
protected:
	int cursor = 0;

public:
	virtual double Sample() = 0;

	virtual double operator()() {
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

	double Sample() override {
		//assert_info(sampler != nullptr, "null sampler");
		double ret = sampler->Sample(cursor++, instance);
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

	double Sample(int d, long long i) {
		return uniform(rng);
	}
private:
	std::mt19937_64 rng;
	std::uniform_real_distribution<double> uniform;
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


class HaltonSampler : public Sampler {
public:
	double Sample(int d, long long i) {
		ASSERT(d < primeList.GetPrimesNum());
		double val = hal(d, i + 1);  // The first one is evil...
		return val;
	}

private:
	inline int rev(const int i, const int p) const {
		return i == 0 ? i : p - i;
	}

	double hal(const int d, long long j) const {
		const int p = primeList.GetPrime(d);
		double h = 0.0, f = 1.0 / p, fct = f;
		while (j > 0) {
			h += rev(j % p, p) * fct;
			j /= p;
			fct *= f;
		}
		return h;
	}
	static PrimeList primeList;
};
PrimeList HaltonSampler::primeList;


struct Vec {
	double x, y, z; // vector: position, also color (r,g,b)
	Vec(double x_ = 0, double y_ = 0, double z_ = 0) { x = x_; y = y_; z = z_; }
	inline Vec operator+(const Vec &b) const { return Vec(x + b.x, y + b.y, z + b.z); }
	inline Vec operator-(const Vec &b) const { return Vec(x - b.x, y - b.y, z - b.z); }
	inline Vec operator+(double b) const { return Vec(x + b, y + b, z + b); }
	inline Vec operator-(double b) const { return Vec(x - b, y - b, z - b); }
	inline Vec operator*(double b) const { return Vec(x * b, y * b, z * b); }
	inline Vec operator*(const Vec &b) const { return Vec(x * b.x, y * b.y, z * b.z); }
	inline Vec operator/(double b) const { if (b == 0) return Vec(); else return Vec(x / b, y / b, z / b); }
	inline bool operator==(const Vec &b) const { return x == b.x && y == b.y && z == b.z; }
	inline bool operator!=(const Vec &b) const { return x != b.x || y != b.y || z != b.z; }
	inline Vec mul(const Vec &b) const { return Vec(x * b.x, y * b.y, z * b.z); }
	inline Vec norm() { return (*this) * (1.0 / sqrt(x * x + y * y + z * z)); }
	inline double dot(const Vec &b) const { return x * b.x + y * b.y + z * b.z; }
	inline double length() const {return sqrt(x * x + y * y + z * z);}
	Vec operator%(Vec&b) const { return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }
	double& operator[](int i) { return i == 0 ? x : i == 1 ? y : z; }
	double maxValue() const {
		return std::max(x, std::max(y, z));
	}
	double Y() const {
		const double YWeight[3] = { 0.212671f, 0.715160f, 0.072169f };
		return YWeight[0] * x + YWeight[1] * y + YWeight[2] * z;
	}
};

Vec operator*(double a, Vec b) { return Vec(a * b.x, a * b.y, a * b.z); }

std::ostream& operator<<(std::ostream &os, const Vec &v) {
	os << v.x << " " << v.y << " " << v.z;
	return os;
}
struct Intersection {
	Vec hit, n, nl, wo;
};

struct AABB {
	Vec minPoint, maxPoint; // axis aligned bounding box
	inline void fit(const Vec &p)
	{
		if (p.x < minPoint.x) minPoint.x = p.x; // min
		if (p.y < minPoint.y) minPoint.y = p.y; // min
		if (p.z < minPoint.z) minPoint.z = p.z; // min
		maxPoint.x = std::max(p.x, maxPoint.x);
		maxPoint.y = std::max(p.y, maxPoint.y);
		maxPoint.z = std::max(p.z, maxPoint.z);
	}
	inline void reset() {
		minPoint = Vec(1e20, 1e20, 1e20);
		maxPoint = Vec(-1e20, -1e20, -1e20);
	}
};

struct HPoint {
	Vec importance, pos, nrm, flux, outDir;
	double r2;
	long long n; // n = N / ALPHA in the paper
	long long pix;
	bool used;
};

struct Ray { Vec o, d; Ray() {}; Ray(Vec o_, Vec d_) : o(o_), d(d_) {} };

void CoordinateSystem(const Vec &v1, Vec *v2, Vec *v3) {
	if (std::abs(v1.x) > std::abs(v1.y))
		*v2 = Vec(-v1.z, 0, v1.x) / std::sqrt(v1.x * v1.x + v1.z * v1.z);
	else
		*v2 = Vec(0, v1.z, -v1.y) / std::sqrt(v1.y * v1.y + v1.z * v1.z);
	*v3 = v1 % (*v2);
}


Vec ConcentricSampleDisk(const Vec &u) {
	// Map uniform Random numbers to $[-1,1]^2$
	Vec uOffset = 2.f * u - Vec(1, 1);

	// Handle degeneracy at the origin
	if (uOffset.x == 0 && uOffset.y == 0) return Vec(0, 0);

	// Apply concentric mapping to point
	double theta, r;
	if (std::abs(uOffset.x) > std::abs(uOffset.y)) {
		r = uOffset.x;
		theta = PiOver4 * (uOffset.y / uOffset.x);
	}
	else {
		r = uOffset.y;
		theta = PiOver2 - PiOver4 * (uOffset.x / uOffset.y);
	}
	return r * Vec(std::cos(theta), std::sin(theta));
}

Vec UniformSampleSphere(const Vec &u) {
	double z = 1 - 2 * u.x;
	double r = std::sqrt(std::max((double)0, (double)1 - z * z));
	double phi = 2 * PI * u.y;
	return Vec(r * std::cos(phi), r * std::sin(phi), z);
}

Vec CosineSampleHemisphere(const Vec &u) {
	Vec d = ConcentricSampleDisk(u);
	double z = std::sqrt(std::max((double)0, 1 - d.x * d.x - d.y * d.y));
	return Vec(d.x, d.y, z);
}

double CosineHemispherePdf(double cosTheta) { return cosTheta * INV_PI; }

class BSDF {
public:
	BSDF(const Intersection &isect) : n(isect.n), nl(isect.nl) {}
	virtual double Pdf(const Vec &wo, const Vec &wi) const = 0;
	virtual Vec Sample_f(const Vec &wo, Vec *wi, double *pdf, Vec rand = Vec(0, 0, 0)) const = 0;
	virtual Vec f(const Vec &wo, const Vec &wi) const { return Vec(0, 0, 0); }
	virtual bool IsDelta() const { return false; }
protected:
	const Vec n, nl;
};

class DiffuseBSDF : public BSDF {
public:
	DiffuseBSDF(const Intersection &isect, Vec r) : BSDF(isect), R(r) {}

	double Pdf(const Vec &wo, const Vec &wi) const {
		return std::abs(wi.dot(nl)) * INV_PI;
	}

	Vec Sample_f(const Vec &wo, Vec *wi, double *pdf, Vec rand) const {
		double r1 = 2. * PI * rand[0], r2 = rand[1];
		double r2s = sqrt(r2);
		Vec w = nl, u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).norm();
		Vec v = w % u;
		*wi = (u* cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();
		*pdf = Pdf(wo, *wi);
		return f(wo, *wi);
	}

	Vec f(const Vec &wo, const Vec &wi) const {
		return R * INV_PI;
	}
private:
	Vec R;
};

class SpecularBSDF : public BSDF {
public:
	SpecularBSDF(const Intersection &isect, Vec r = Vec(1.0, 1.0, 1.0)) : BSDF(isect), R(r) {}

	double Pdf(const Vec &wo, const Vec &wi) const {
		return 0.0;
	}

	Vec Sample_f(const Vec &wo, Vec *wi, double *pdf, Vec rand) const {
		*wi = (nl * 2.0 * nl.dot(wo) - wo).norm();
		*pdf = 1.0;
		double cosTheta = std::abs((*wi).dot(n));
		return R / cosTheta;
	}

	Vec f(const Vec &wo, const Vec &wi) const {
		return Vec();
	}

	bool IsDelta() const { return true; }
private:
	Vec R;
};

class TransmissionBSDF : public BSDF {
public:
	TransmissionBSDF(const Intersection &isect, Vec fa = Vec(1.0, 1.0, 1.0), double eta1 = 1.0, double eta2 = 1.5) :
		BSDF(isect), Fa(fa), nc(eta1), nt(eta2) {}

	double Pdf(const Vec &wo, const Vec &wi) const {
		return 0.0;
	}

	Vec Sample_f(const Vec &wo, Vec *wi, double *pdf, Vec rand) const {
		bool into = (n.dot(nl) > 0.0);
		double nnt = into ? nc / nt : nt / nc, ddn = (-1 * wo).dot(nl), cos2t;
		// total internal reflection
		if ((cos2t = 1 - nnt * nnt*(1 - ddn * ddn)) < 0) {
			*wi = (nl * 2.0 * nl.dot(wo) - wo).norm();
			double cosTheta = std::abs((*wi).dot(n));
			*pdf = 1.0;
			//return Fa / cosTheta;
			return Vec();
		}
		Vec td = ((-1 * wo) * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
		double Re = Fresnell(wo, td, n, nl);
		double P = Re * 0.5 + 0.25;
		if (rand.z < P) {
			*wi = (nl * 2.0 * nl.dot(wo) - wo).norm();
			*pdf = P;
			double cosTheta = std::abs((*wi).dot(n));
			return Fa * Re / cosTheta;
		}
		else {

			*wi = td;
			*pdf = 1 - P;
			double cosTheta = std::abs((*wi).dot(n));
			return Fa * (1.0 - Re) / cosTheta;
		}
	}

	Vec f(const Vec &wo, const Vec &wi) const {
		return Vec();
	}

	double Fresnell(const Vec &wo, const Vec &td, const Vec &n, const Vec &nl) const {
		bool into = (n.dot(nl) > 0.0);
		double nnt = into ? nc / nt : nt / nc, ddn = (-1 * wo).dot(nl);//, cos2t;
		//cos2t = 1 - nnt * nnt*(1 - ddn * ddn);
		double a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : td.dot(n));
		double Re = R0 + (1 - R0) * c * c* c * c * c;
		return Re;
	}

	double FrDielectric(double cosThetaI, double etaI, double etaT) {
		cosThetaI = Clamp(cosThetaI, -1, 1);
		// Potentially swap indices of refraction
		bool entering = cosThetaI > 0.f;
		if (!entering) {
			std::swap(etaI, etaT);
			cosThetaI = std::abs(cosThetaI);
		}

		// Compute _cosThetaT_ using Snell's law
		double sinThetaI = std::sqrt(std::max((double)0, 1 - cosThetaI * cosThetaI));
		double sinThetaT = etaI / etaT * sinThetaI;

		// Handle total internal reflection
		if (sinThetaT >= 1) return 1;
		double cosThetaT = std::sqrt(std::max((double)0, 1 - sinThetaT * sinThetaT));
		double Rparl = ((etaT * cosThetaI) - (etaI * cosThetaT)) /
			((etaT * cosThetaI) + (etaI * cosThetaT));
		double Rperp = ((etaI * cosThetaI) - (etaT * cosThetaT)) /
			((etaI * cosThetaI) + (etaT * cosThetaT));
		return (Rparl * Rparl + Rperp * Rperp) / 2;
	}

	bool IsDelta() const { return true; }

private:
	double nc, nt;
	Vec Fa;
};

class Shape {
public:
	Shape(ReflectionType type, const Vec &color, const Vec &emission, bool isL = false): 
		reflType(type), c(color), e(emission), isLight(isL){}
	virtual double Intersect(const Ray &r, Intersection *isect) const = 0;
	virtual Vec Sample(double *pdf, Vec rand) const = 0;
	virtual Vec Sample(const Intersection &isect, double *pdf, Vec u) const = 0;
	virtual std::shared_ptr<BSDF> GetBSDF(const Intersection &isect) const {
		if (reflType == DIFF) {
			return std::dynamic_pointer_cast<BSDF>(std::make_shared<DiffuseBSDF>(isect, c));
		}
		else if (reflType == SPEC) {
			return std::dynamic_pointer_cast<BSDF>(std::make_shared<SpecularBSDF>(isect, c));
		}
		else if (reflType == REFR) {
			return std::dynamic_pointer_cast<BSDF>(std::make_shared<TransmissionBSDF>(isect, c));
		}
		return std::shared_ptr<BSDF>();
	}
	bool IsLight() const { return isLight; }

	virtual Vec GetNorm(const Vec &point) const = 0;

	virtual Vec GetEmission() const { return e; }

	int GetId() const { return shapeId; }

	virtual double Area() const = 0;

	friend class Scene;
private:
	ReflectionType reflType;
	Vec c, e;
	bool isLight;
	int shapeId;
};

class Light {
public:
	Light() {}
	virtual Vec DirectIllumination(const Intersection &isect, const std::shared_ptr<BSDF> &bsdf, const Vec &importance, Vec *dir, Vec u) const = 0;
	virtual Vec Emission() const = 0;
	virtual Vec SampleLight(Vec *pos, Vec *dir, Vec *lightNorm, double *pdfPos, double *pdfDir, Vec u, Vec v) const {
		SampleOnLight(pos, dir, lightNorm, pdfPos, pdfDir, u, v);
		return Emission();
	}
	virtual int GetId() const = 0;
	virtual Vec Power() const = 0;
	virtual bool IsAreaLight() const { return false; }
	virtual std::shared_ptr<Shape> GetShapePtr() const = 0;
protected:
	virtual void SampleOnLight(Vec *pos, Vec *dir, Vec *lightNorm, double *pdfPos, double *pdfDir, Vec u, Vec v) const = 0;
};


class Sphere: public Shape {
public:
	Sphere(double radius, Vec position, Vec emission, Vec color, ReflectionType reflType) : 
		rad(radius), p(position), Shape(reflType, color, emission, emission != Vec()){}

	double Intersect(const Ray &r, Intersection *isect) const {
		// ray-sphere Intersection returns distance
		Vec op = p - r.o;
		double t, b = op.dot(r.d), det = b * b - op.dot(op) + rad * rad;
		if (det < 0) {
			return 1e20;
		}
		else {
			det = sqrt(det);
		}
		t = (t = b - det) > 1e-4 ? t : ((t = b + det) > 1e-4 ? t : 1e20);

		isect->hit = r.o + r.d * t;
		isect->n = (isect->hit - p).norm();
		isect->nl = isect->n.dot(r.d) < 0 ? isect->n : isect->n * -1;
		isect->wo = -1 * r.d;
		return t;
	}

	Vec Sample(double *pdf, Vec rand) const {
		*pdf = 1.0 / (4.0 * PI * rad * rad);
		return UniformSampleSphere(rand) * rad + p;
	}

	Vec Sample(const Intersection &isect, double *pdf, Vec u) const {
		Vec sw = p - isect.hit, su = ((fabs(sw.x) > .1 ? Vec(0, 1) : Vec(1)) % sw).norm(), sv = sw % su;
		double cos_a_max = sqrt(1 - rad * rad / (isect.hit - p).dot(isect.hit - p));
		double zeta1 = u.x, zeta2 = u.y;
		double cos_a = 1 - zeta1 + zeta1 * cos_a_max;
		double sin_a = sqrt(1 - cos_a * cos_a);
		double phi = 2 * PI * zeta2;
		Vec dir = su * cos(phi) * sin_a + sv * sin(phi) * sin_a + sw * cos_a;
		double omega = 2 * PI *(1 - cos_a_max);
		*pdf = 1.0 / omega;
		return dir.norm();
	}

	Vec GetNorm(const Vec & point) const {
		return (point - p).norm();
	}

	double Area() const {
		return 4.0 * PI * rad * rad;
	}

private:
	double rad; Vec p;
};


class AreaLight : public Light{
public:
	AreaLight(const std::shared_ptr<Shape>& pShape): shape(pShape){}
	Vec DirectIllumination(const Intersection &isect, const std::shared_ptr<BSDF> &bsdf, const Vec &importance, Vec *dir, Vec u) const {
		double pdf;
		*dir = shape->Sample(isect, &pdf, u);
		Vec f = bsdf->f(isect.wo, *dir);
		return importance * f * std::abs((*dir).dot(isect.nl)) * Emission() / pdf; 
	}

	Vec Emission() const {
		return shape->GetEmission();
	}

	int GetId() const {
		return shape->GetId();
	}

	Vec Power() const {
		return Emission() * shape->Area() * PI;
	}

	bool IsAreaLight() const { return true; }

	std::shared_ptr<Shape> GetShapePtr() const { return shape; }
protected:
	void SampleOnLight(Vec *pos, Vec *dir, Vec *lightNorm, double *pdfPos, double *pdfDir, Vec u, Vec v) const {
		//sample a position
		*pos = shape->Sample(pdfPos, u);
		*lightNorm = shape->GetNorm(*pos);
		Vec ss, ts;
		CoordinateSystem(*lightNorm, &ss, &ts);
		Vec dirLocal = CosineSampleHemisphere(v);
		double cosTheta = dirLocal.z;
		*dir = (ss * dirLocal.x + ts * dirLocal.y + *lightNorm * dirLocal.z).norm();
		*pdfDir = CosineHemispherePdf(cosTheta);
	}
private:
	std::shared_ptr<Shape> shape;
};

class Filter
{
protected:
	double radius;

public:
	Filter(const double rad)
		: radius(rad)
	{
	}
	virtual ~Filter() {}

	const double GetRadius() const
	{
		return radius;
	}
	virtual double Evaluate(const double dx, const double dy) const = 0;
};

class BoxFilter : public Filter
{
public:
	BoxFilter()
		: Filter(0.5)
	{
	}

public:
	double Evaluate(const double dx, const double dy) const {
		return 1.0;
	}
};

class Film {
protected:
	struct Pixel
	{
		Vec color;
		Vec splat;
		double weight;
	};
public:
	Film(int w, int h) : resX(w), resY(h) {
		aspect = (double)(resX) / (double)(resY);
		filter = std::unique_ptr<BoxFilter>(new BoxFilter());
	}
	double Area() const {
		return area;
	}
	void AddSample(double x, double y, const Vec &sample) {
		x -= 0.5;
		y -= 0.5;

		int minX = (int)(std::ceil(x - filter->GetRadius()));
		int maxX = (int)(std::floor(x + filter->GetRadius()));
		int minY = (int)(std::ceil(y - filter->GetRadius()));
		int maxY = (int)(std::floor(y + filter->GetRadius()));
		minX = std::max(0, minX);
		maxX = std::min(maxX, resX - 1);
		minY = std::max(0, minY);
		maxY = std::min(maxY, resY - 1);

		for (int i = minY; i <= maxY; i++) {
			for (int j = minX; j <= maxX; j++) {
				//int rowAdd = resY - 1 - i;
				//int colAdd = j;
				int pixelIndex = i * resX + j;
				Pixel& pixel = pixelBuffer[pixelIndex];

				double weight = filter->Evaluate(j - x, i - y);
				pixel.weight += weight;
				pixel.color = pixel.color + sample * weight;
			}
		}
	}

	void Film::AddSplat(double x, double y, const Vec& sample)
	{
		int X = (int)(std::floor(x));
		int Y = (int)(std::floor(y));
		X = Clamp(X, 0, resX - 1);
		Y = Clamp(Y, 0, resY - 1);

		//int rowAdd = resY - 1 - Y;
		//int colAdd = X;
		int pixelIndex = X * resX + Y;
		Pixel& pixel = pixelBuffer[pixelIndex];
		pixel.splat = pixel.splat + sample;
	}

	void SetImage(const std::vector<Vec> &image) {
		imageBuffer = image;
	}

	void SetFileName(const std::string &pFileName) {
		filename = pFileName;
	}

	void SaveImage() {
		std::string suffix = filename.substr(filename.size() - 4);
		if (suffix == ".png") {
			WritePngFile();
		}
		else if (suffix == ".bmp") {
			WriteBmpFile();
		}
		else {

		}
	}
public:
	int resX, resY;
	double width, heigh;
	double aspect;
	double area;
	Vec LL, LU, RL, RU;
private:
	struct BmpHeader
	{
		unsigned int   mFileSize;        // Size of file in bytes
		unsigned int   mReserved01;      // 2x 2 reserved bytes
		unsigned int   mDataOffset;      // Offset in bytes where data can be found (54)

		unsigned int    mHeaderSize;      // 40B
		unsigned int    mWidth;           // Width in pixels
		unsigned int    mHeight;          // Height in pixels

		short  mColorPlates;     // Must be 1
		short  mBitsPerPixel;    // We use 24bpp
		unsigned int   mCompression;     // We use BI_RGB ~ 0, uncompressed
		unsigned int   mImageSize;       // mWidth x mHeight x 3B
		unsigned int   mHorizRes;        // Pixels per meter (75dpi ~ 2953ppm)
		unsigned int   mVertRes;         // Pixels per meter (75dpi ~ 2953ppm)
		unsigned int   mPaletteColors;   // Not using palette - 0
		unsigned int   mImportantColors; // 0 - all are important
	};

	void WritePngFile() {
		FILE* f = fopen(filename.c_str(), "wb");
		typedef unsigned char byte;
		byte *pngImage = new byte[resX * resY * 3];
		for (int j = 0; j < resY; ++j) {
			for (int i = 0; i < resX; ++i) {
				int index = 3 * j * resX + 3 * i;
				int imageBufferIndex = j * resX + i;
				pngImage[index] = (byte)(toInt(imageBuffer[imageBufferIndex].x));
				pngImage[index + 1] = (byte)(toInt(imageBuffer[imageBufferIndex].y));
				pngImage[index + 2] = (byte)(toInt(imageBuffer[imageBufferIndex].z));
			}
		}
		svpng(f, resX, resY, pngImage, 0);
		delete[] pngImage;
		fclose(f);
	}

	// tone mapping and gamma correction
	//int toInt(double x) {
	//	return int(pow(1 - exp(-x), 1 / 2.2) * 255 + .5);
	//}
	int toInt(double x) { 
		return int(pow(Clamp(x, 0.0, 1.0), 1 / 2.2) * 255 + .5); 
	}

	void WriteToPixelBuffer() {
		for (int i = 0; i < resX; ++i) {
			for (int j = 0; j < resY; ++j) {
				int index = i * resX + j;
				Pixel &pixel = pixelBuffer[index];
				imageBuffer[index] = pixel.color / pixel.weight + pixel.splat;
			}
		}
	}

	void WriteBmpFile()  
	{
		std::ofstream bmp(filename.c_str(), std::ios::binary);
		BmpHeader header;
		bmp.write("BM", 2);
		header.mFileSize = (unsigned int)(sizeof(BmpHeader) + 2) + resX * resY * 3;
		header.mReserved01 = 0;
		header.mDataOffset = (unsigned int)(sizeof(BmpHeader) + 2);
		header.mHeaderSize = 40;
		header.mWidth = resX;
		header.mHeight = resY;
		header.mColorPlates = 1;
		header.mBitsPerPixel = 24;
		header.mCompression = 0;
		header.mImageSize = resX * resY * 3;
		header.mHorizRes = 2953;
		header.mVertRes = 2953;
		header.mPaletteColors = 0;
		header.mImportantColors = 0;
		bmp.write((char*)&header, sizeof(header));

		for (int y = 0; y < resY; y++)
		{
			for (int x = 0; x < resX; x++)
			{
				const Vec &rgbF = imageBuffer[x + (resY - y - 1) * resX];
				typedef unsigned char byte;
				byte bgrB[3];
				bgrB[0] = byte(toInt(rgbF.z));
				bgrB[1] = byte(toInt(rgbF.y));
				bgrB[2] = byte(toInt(rgbF.x));

				bmp.write((char*)&bgrB, sizeof(bgrB));
			}
		}
	}

private:
	std::string filename;
	std::vector<Pixel> pixelBuffer;
	std::vector<Vec> imageBuffer;
	std::unique_ptr<Filter> filter;
	std::vector<Spinlock> bufferLocks;

};

class Camera {
public:
	Camera(const std::shared_ptr<Film> &pFilm) {
		film = pFilm;
	}
	virtual Ray GenerateRay(int pixelX, int pixelY, Vec u, double offset = 0.0) const = 0;
	virtual Vec We(const Ray &ray) const = 0;
	virtual Vec Sample_Wi(const Intersection &isect, double *pdfW, Vec *wi, Vec u) const = 0;
	virtual double PdfPos() const = 0;
	virtual double PdfDir(const Ray &cameraRay) const = 0;
	virtual std::shared_ptr<Film> GetFilm() const { return film; }
protected:
	std::shared_ptr<Film> film;
};

class PinHoleCamera : public Camera {
public:
	PinHoleCamera(const std::shared_ptr<Film> &pFilm, const Vec &position, 
		const Vec &pCz, const Vec &pCx, const Vec &pCy, double pFovy, double dis) :
		Camera(pFilm), pos(position), cz(pCz), cx(pCx), cy(pCy), fovy(pFovy), filmDistance(dis) {
		Initialize();
	}

	void Initialize() {
		Vec filmCenter = pos + cz * filmDistance;
		//std::cout << "film cetner: " << filmCenter << std::endl;
		film->heigh = filmDistance * std::tan(fovy * 0.5 * PI / 180) * 2.0;
		film->width = film->heigh * film->aspect;
		film->area = film->width * film->heigh;
		/*
		film->LU = filmCenter + cy * film->heigh * 0.5 - cx * film->width * 0.5;
		film->LL = filmCenter - cy * film->heigh * 0.5 - cx * film->width * 0.5;
		film->RU = filmCenter + cy * film->heigh * 0.5 + cx * film->width * 0.5;
		film->RL = filmCenter - cy * film->heigh * 0.5 + cx * film->width * 0.5;
		*/
	}

	Ray GenerateRay(int pixelX, int pixelY, Vec u, double offset) const {
		
		Vec dir = cx * ((pixelX + u.x) / film->resX - 0.5) * film->width + 
			cy * (-(pixelY + u.y) / film->resY + 0.5) * film->heigh  + cz * filmDistance;
		dir = dir.norm();

		return Ray(pos + dir * offset, dir);
	}

	Vec We(const Ray &ray) const {
		double pdfA = 1.0; // for the pinhole camera
		double area = film->area;
		double cosTheta = cz.dot(ray.d);
		double cos2Theta = cosTheta * cosTheta;
		double value = filmDistance * filmDistance * pdfA / (area * cos2Theta * cos2Theta);
		return Vec(value, value, value);
	}

	Vec Sample_Wi(const Intersection &isect, double *pdfW, Vec *wi, Vec u) const {
		*wi = (pos - isect.hit);
		double distance = wi->length();
		*wi = wi->norm();
		double cosTheta = cz.dot(-1 * (*wi));
		*pdfW = 1.0 * (distance * distance) / cosTheta;
		//*PdfW = 1.0 * (dis / CosTheta) * (dis / CosTheta) / CosTheta;
		return We(Ray(isect.hit, -1 * (*wi)));
	}

	double PdfPos() const {
		return 1.0;
	}

	double PdfDir(const Ray &cameraRay) const {
		double filmArea = film->Area();
		double cosTheta = std::abs(cz.dot(cameraRay.d));
		double cos2Theta = cosTheta * cosTheta;
		return filmDistance * filmDistance / (filmArea * cos2Theta * cosTheta);
	}
private:
	Vec pos, cx, cy, cz;
	double fovy;
	double filmDistance;

};

//std::vector<double> radius2;
//std::vector<long long> photonNums;
//std::vector<Vec> flux;
//std::vector<Vec> directillum;
//std::vector<HPoint> hitPoints;
//std::vector<std::vector<HPoint*>> hashGrid;
//std::vector<Spinlock> hashGridSpinlocks;
//long long hashNum, pixelIndex;
//double hashCellSize;
//AABB hpbbox;



class HashGrid {
public:
	HashGrid(double initialRadius = 1.0) {
		irad = initialRadius;
	}

	void Initialize(double initialRadius) {
		irad = initialRadius;
	}

	// spatial hash function
	inline unsigned int hash(const int ix, const int iy, const int iz) {
		return (unsigned int)((ix * 73856093) ^ (iy * 19349663) ^ (iz * 83492791)) % hashNum;
	}

	void BuildHashGrid(std::vector<HPoint> &hitPoints) {
		// find the bounding box of all the measurement points
		hpbbox.reset();
		for (const HPoint &hp : hitPoints) {
			if (hp.used) hpbbox.fit(hp.pos);
		}

		// heuristic for initial radius
		Vec ssize = hpbbox.maxPoint - hpbbox.minPoint;
		//double irad = ((ssize.x + ssize.y + ssize.z) / 3.0) / ((w + h) / 2.0) * 2.0;
		double irad = 0.8;
		// determine hash table size
		// we now find the bounding box of all the measurement points inflated by the initial radius
		hpbbox.reset();
		int vphoton = 0;

		for (const HPoint &hp : hitPoints) {
			if (!hp.used) continue;
			vphoton++;
			hpbbox.fit(hp.pos - irad);
			hpbbox.fit(hp.pos + irad);
		}

		// make each grid cell two times larger than the initial radius
		invHashCellSize = 1.0 / (irad * 2.0);
		hashNum = vphoton;

		// build the hash table

		hashGrid.resize(hashNum);
		hashGridSpinlocks.resize(hashNum);


		//for (HPoint &hp : hitPoints) {
		int hitPointNum = (int)hitPoints.size();
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < hitPointNum; ++i) {
			HPoint &hp = hitPoints[i];
			if (!hp.used) continue;
			Vec BMin = ((hp.pos - irad) - hpbbox.minPoint) * invHashCellSize;
			Vec BMax = ((hp.pos + irad) - hpbbox.minPoint) * invHashCellSize;
			for (int iz = abs(int(BMin.z)); iz <= abs(int(BMax.z)); iz++)
			{
				for (int iy = abs(int(BMin.y)); iy <= abs(int(BMax.y)); iy++)
				{
					for (int ix = abs(int(BMin.x)); ix <= abs(int(BMax.x)); ix++)
					{
						int hv = hash(ix, iy, iz);
						std::lock_guard<Spinlock> lock(hashGridSpinlocks[hv]);
						hashGrid[hv].push_back(&hp);
					}
				}
			}
		}
	}

	std::vector<HPoint*>& GetGrid(const Vec &point) {
		// photon ray
		// find neighboring measurement points and accumulate flux via progressive density estimation
		Vec hh = (point - hpbbox.minPoint) * invHashCellSize;
		int ix = std::abs(int(hh.x)), iy = std::abs(int(hh.y)), iz = std::abs(int(hh.z));
		// strictly speaking, we should use #pragma omp critical here.
		// it usually works without an artifact due to the fact that photons are 
		// rarely accumulated to the same measurement points at the same time (especially with QMC).
		// it is also significantly faster.
		return hashGrid[hash(ix, iy, iz)];
	}

	void ClearHashGrid() {
		for (auto &e : hashGrid) {
			std::vector<HPoint*>().swap(e);
		}
	}

private:
	std::vector<std::vector<HPoint*>> hashGrid;
	std::vector<Spinlock> hashGridSpinlocks;
	long long hashNum, pixelIndex;
	double invHashCellSize;
	double irad;
	AABB hpbbox;
};

/*
Sphere sph[] = { // Scene: radius, position, color, material
	Sphere(1e5, Vec(1e5 + 1,40.8,81.6), Vec(), Vec(.75,.25,.25),DIFF),//Left
	Sphere(1e5, Vec(-1e5 + 99,40.8,81.6), Vec(), Vec(.25,.25,.75),DIFF),//Right
	Sphere(1e5, Vec(50,40.8, 1e5), Vec(),Vec(.75,.75,.75),DIFF),//Back
	Sphere(1e5, Vec(50,40.8,-1e5 + 170), Vec(), Vec(), DIFF),//Front
	Sphere(1e5, Vec(50, 1e5, 81.6),  Vec(), Vec(.75,.75,.75),DIFF),//Bottomm
	Sphere(1e5, Vec(50,-1e5 + 81.6,81.6),Vec(),Vec(.75,.75,.75),DIFF),//Top
	Sphere(16.5,Vec(27,16.5,47), Vec(), Vec(1,1,1)*.999, SPEC),//Mirror
	Sphere(16.5,Vec(73,16.5,88), Vec(), Vec(1,1,1)*.999, REFR),//Glass
	Sphere(8.5, Vec(50,8.5,60),  Vec(), Vec(1,1,1)*.999, DIFF),
	Sphere(6.5, Vec(73,16.5,88),  Vec(), Vec(1,1,1)*.999, DIFF),
	Sphere(8.0, Vec(50,81.6 - 16.5,81.6),Vec(0.25,0.25,0.25) * 100,  Vec(), DIFF),//Lite
};//Middle*/
/*
Sphere sph[] = {//Scene: radius, position, emission, color, material
	Sphere(1e5, Vec(1e5 + 1,40.8,81.6), Vec(),Vec(.75,.25,.25),DIFF),//Left
	Sphere(1e5, Vec(-1e5 + 99,40.8,81.6),Vec(),Vec(.25,.25,.75),DIFF),//Rght
	Sphere(1e5, Vec(50,40.8, 1e5),     Vec(),Vec(.75,.75,.75),DIFF),//Back
	Sphere(1e5, Vec(50,40.8,-1e5 + 170), Vec(),Vec(.75,.75,.75), DIFF),//Frnt
	Sphere(1e5, Vec(50, 1e5, 81.6),    Vec(),Vec(.75,.75,.75),DIFF),//Botm
	Sphere(1e5, Vec(50,-1e5 + 81.6,81.6),Vec(),Vec(.75,.75,.75),DIFF),//Top
	Sphere(16.5,Vec(27,16.5,47),       Vec(),Vec(1,1,1)*.999, SPEC),//Mirr
	Sphere(16.5,Vec(73,16.5,78),       Vec(),Vec(1,1,1)*.999, REFR),//Glas
	Sphere(8.0, Vec(50,81.6 - 16.5,81.6),Vec(0.3,0.3,0.3) * 100,  Vec(), DIFF),//Lite
};*/

Sphere sph[] = {//Scene: radius, position, emission, color, material
	Sphere(1e5, Vec(1e5 + 1,40.8,81.6), Vec(),Vec(.75,.25,.25),DIFF),//Left
	Sphere(1e5, Vec(-1e5 + 99,40.8,81.6),Vec(),Vec(.25,.25,.75),DIFF),//Rght
	Sphere(1e5, Vec(50,40.8, 1e5),     Vec(),Vec(.75,.75,.75),DIFF),//Back
	Sphere(1e5, Vec(50,40.8,-1e5 + 170), Vec(),Vec(.75,.75,.75), DIFF),//Frnt
	Sphere(1e5, Vec(50, 1e5, 81.6),    Vec(),Vec(.75,.75,.75),DIFF),//Botm
	Sphere(1e5, Vec(50,-1e5 + 81.6,81.6),Vec(),Vec(.75,.75,.75),DIFF),//Top
	Sphere(16.5,Vec(27,16.5,47),       Vec(),Vec(1,1,1)*.999, REFR),//Mirr
	Sphere(7.0,Vec(27,16.5,47),       Vec(),Vec(.25,.25,.75), DIFF),//Mirr
	Sphere(16.5,Vec(73,26.5,78),       Vec(),Vec(1,1,1)*.999, REFR),//Glas
	Sphere(9.5, Vec(53,9.5,88),  Vec(), Vec(1,1,1)*.999, REFR),
	Sphere(8.0, Vec(50,81.6 - 16.5,81.6),Vec(0.3,0.3,0.3) * 100,  Vec(), DIFF),//Lite
};


class Scene {
public:

	void SetCamera(const std::shared_ptr<Camera> &pCamera) {
		camera = pCamera;
		shapeNum = 0;
	}

	void AddShape(std::shared_ptr<Shape> shape) {
		shapes.push_back(shape);
		shapes[shapeNum]->shapeId = shapeNum;
		++shapeNum;
	}

	void AddLight(std::shared_ptr<Light> light) {
		lights.push_back(light);
		if (light->IsAreaLight()) {
			AddShape(light->GetShapePtr());
		}
	}

	void Initialize() {
		lightPowerDistribution = ComputeLightPowerDistribution();
	}

	bool Intersect(const Ray &r, double *t, Intersection *isect, std::shared_ptr<Shape> &hitObj) const {
		int n = (int)shapes.size();
		double d;
		*t = Inf;
		for (int i = 0; i < n; ++i) {
			Intersection intersection;
			d = shapes[i]->Intersect(r, &intersection);
			if (d < *t) {
				*t = d;
				hitObj = shapes[i];
				*isect = intersection;
			}
		}
		return *t < Inf;
	}

	const std::vector<std::shared_ptr<Light>>& GetLights() const {
		return lights;
	}

	const std::vector<std::shared_ptr<Shape>>& GetShapes() const {
		return shapes;
	}

	std::shared_ptr<Camera> GetCamera() const {
		return camera;
	}

	std::shared_ptr<Light> SampleOneLight(double *lightPdf, double u) const {
		int nLights = (int)(lights.size());
		int lightNum;
		if (lightPowerDistribution != nullptr) {
			lightNum = lightPowerDistribution->SampleDiscrete(u, lightPdf);
			if (*lightPdf == 0) return nullptr;
		} else {
			lightNum = std::min((int)(u * nLights), nLights - 1);
			*lightPdf = 1.0 / nLights;
		}
		return lights[lightNum];
	}

private:

	std::unique_ptr<Distribution1D> ComputeLightPowerDistribution() {
		if (lights.size() == 0) return nullptr;
		std::vector<double> lightPower;
		for (const auto &light : lights)
			lightPower.push_back(light->Power().Y());
		return std::unique_ptr<Distribution1D>(
			new Distribution1D(&lightPower[0], (int)lightPower.size()));
	}

private:
	std::vector<std::shared_ptr<Shape>> shapes;
	std::vector<std::shared_ptr<Light>> lights;
	std::shared_ptr<Camera> camera;
	std::unique_ptr<Distribution1D> lightPowerDistribution;
	int shapeNum;
};


class Integrator {
public:
	virtual void Render(const Scene &scene) = 0;
	virtual ~Integrator() {}
	static Vec DirectIllumination(const Scene &scene, const Intersection &isect, const std::shared_ptr<BSDF> &bsdf, Vec importance, Vec u) {
		Vec L;
		const std::vector<std::shared_ptr<Light>> lights = scene.GetLights();
		for (auto light : lights) {
			Vec dir;
			std::shared_ptr<Shape> hitObj;
			double t;
			Vec Li = light->DirectIllumination(isect, bsdf, importance, &dir, u);
			Intersection intersection;
			if (scene.Intersect(Ray(isect.hit, dir), &t, &intersection, hitObj) && hitObj->GetId() == light->GetId()) {
				L = L + Li;
			}
		}
		return L;
	}
};


class SPPM : public Integrator {
public:
	SPPM(int iterations, int nPhotonsPerStage, int maxDepth, double initialRadius, const std::shared_ptr<Sampler> &pSampler):
		nIterations(iterations), nPhotonsPerRenderStage(nPhotonsPerStage), 
		maxDepth(maxDepth), initialRadius(initialRadius), alpha(ALPHA), sampler(pSampler)
	{
		
	}

	void GenratePhoton(const Scene &scene, Ray *pr, Vec *f, double u, const Vec &v, const Vec &w) {
		double lightPdf;
		std::shared_ptr<Light> light = scene.SampleOneLight(&lightPdf, u);
		Vec Le = light->Emission();
		Vec pos, lightDir, lightNorm;
		double pdfPos, pdfDir;
		light->SampleLight(&pos, &lightDir, &lightNorm, &pdfPos, &pdfDir, v, w);
		pr->o = pos + lightDir * rayeps;
		pr->d = lightDir;
		double cosTheta = std::abs(lightNorm.dot(lightDir));
		*f = Le * cosTheta / (pdfPos * pdfDir * lightPdf);
	}

	void TraceEyePath(const Scene &scene, StateSequence &rand, const Ray &ray, long long pixel) {
		Ray r = ray;
		bool deltaBoundEvent = false;
		Vec importance(1.0, 1.0, 1.0);
		for (int i = 0; i < maxDepth; ++i) {
			double t;
			Intersection isect;
			std::shared_ptr<Shape> hitObj;
			if (!scene.Intersect(r, &t, &isect, hitObj)) return;
			std::shared_ptr<BSDF> bsdf = hitObj->GetBSDF(isect);
			Vec wi;
			double pdf;
			if (bsdf->IsDelta()) {
				Vec f = bsdf->Sample_f(-1 * r.d, &wi, &pdf, Vec(rand(), rand(), rand()));
				importance = f * std::abs(wi.dot(isect.n)) * importance / pdf;
				r.o = isect.hit + wi * rayeps;
				r.d = wi;
				deltaBoundEvent = true;
			}
			else {
				HPoint &hp = hitPoints[pixel];
				hp.used = true;
				hp.importance = importance;
				hp.pos = isect.hit;
				hp.nrm = isect.n;
				hp.pix = pixel;
				hp.outDir = -1 * r.d;
				//hitPoints[pixel] = hp;
				if ((i == 0 || deltaBoundEvent) && hitObj->IsLight())
					directillum[hp.pix] = directillum[hp.pix] + importance * hitObj->GetEmission();
				else
					directillum[hp.pix] = directillum[hp.pix] +
					DirectIllumination(scene, isect, bsdf, hp.importance, Vec(rand(), rand(), rand()));

				return;
			}
		}
	}

	void TracePhoton(const Scene &scene, StateSequence &rand, const Ray &ray, Vec photonFlux) {
		Ray r = ray;
		for (int i = 0; i < maxDepth; ++i) {
			double t;
			Intersection isect;
			std::shared_ptr<Shape> hitObj;
			if (!scene.Intersect(r, &t, &isect, hitObj)) return;
			std::shared_ptr<BSDF> bsdf = hitObj->GetBSDF(isect);
			Vec wi;
			double pdf;
			Vec f = bsdf->Sample_f(-1 * r.d, &wi, &pdf, Vec(rand(), rand(), rand()));
			Vec estimation = f * std::abs(wi.dot(isect.n)) / pdf;
			if (bsdf->IsDelta()) {
				photonFlux = photonFlux * estimation;
				r.o = isect.hit + wi * rayeps;
				r.d = wi;
			}
			else {
				if (i > 0) {
					std::vector<HPoint*> &hp = hashGrid.GetGrid(isect.hit);
					for (HPoint* hitpoint : hp) {
						Vec v = hitpoint->pos - isect.hit;
						if ((hitpoint->nrm.dot(isect.n) > 1e-3) && (v.dot(v) <= radius2[hitpoint->pix])) {
							// unlike N in the paper, hitpoint->n stores "N / ALPHA" to make it an integer value
							double g = (photonNums[hitpoint->pix] * alpha + alpha) / (photonNums[hitpoint->pix] * alpha + 1.0);
							radius2[hitpoint->pix] = radius2[hitpoint->pix] * g;
							photonNums[hitpoint->pix]++;
							Vec contribution = hitpoint->importance * bsdf->f(hitpoint->outDir, -1 * r.d) * photonFlux;
							flux[hitpoint->pix] = (flux[hitpoint->pix] + contribution) * g;
						}
					}
				}
				double p = estimation.maxValue();
				if (p < 1) {
					if (rand() < p) photonFlux = photonFlux / p;
					else break;
				}
				photonFlux = photonFlux * estimation;
				r.o = isect.hit + wi * rayeps;
				r.d = wi;
			}
		}
	}

	void Render(const Scene &scene) {
		fprintf(stderr, "Rendering ...\n");
		int resX = scene.GetCamera()->GetFilm()->resX;
		int resY = scene.GetCamera()->GetFilm()->resY;
		Initialize(resX, resY);
		for (int iter = 0; iter < nIterations; ++iter) {
			std::shared_ptr<Sampler> randomSampler = std::shared_ptr<Sampler>(new RandomSampler(resX * resY));
#pragma omp parallel for schedule(guided)
			for (int y = 0; y < resY; y++) {
				//fprintf(stderr, "\rHitPointPass %5.2f%%", 100.0*y / (h - 1));
				for (int x = 0; x < resX; x++) {
					int pixel = x + y * resX;
					hitPoints[pixel].used = false;
					RandomStateSequence rand(randomSampler, iter);
					//RandomStateSequence rand(sampler, iter);
					Ray ray = scene.GetCamera()->GenerateRay(x, y, Vec(rand(), rand(), rand()), 140);
					TraceEyePath(scene, rand, ray, pixel);
				}
			}
			
			hashGrid.BuildHashGrid(hitPoints);
			{
				//fprintf(stderr, "\n");

				//fprintf(stderr, "\rPhotonPass %5.2f%%", percentage);
				Ray ray;
				Vec photonFlux;
#pragma omp parallel for schedule(guided)
				for (int j = 0; j < nPhotonsPerRenderStage; j++) {
					RandomStateSequence rand(sampler, iter * nPhotonsPerRenderStage + j);
					GenratePhoton(scene, &ray, &photonFlux, rand(), Vec(rand(), rand(), rand()), Vec(rand(), rand(), rand()));
					TracePhoton(scene, rand, ray, photonFlux);
				}
				//fprintf(stderr, "\n");

			}
			hashGrid.ClearHashGrid();
			double percentage = 100.*(iter + 1) / nIterations;
			fprintf(stderr, "\rIterations: %5.2f%%", percentage);
		}

		// density estimation
		std::cout << "\nflux size: " << flux.size() << std::endl;
#pragma omp parallel for schedule(guided)
		for (int i = 0; i < flux.size(); ++i) {
			c[i] = c[i] + flux[i] * (1.0 / (PI * radius2[i] * nIterations * nPhotonsPerRenderStage))
				+ directillum[i] / nIterations;
			//c[i] = c[i] + directillum[i] / nIterations;
		}
		scene.GetCamera()->GetFilm()->SetImage(c);
	}

private:
	void Initialize(int w, int h) {
		radius2.resize(w * h);
		photonNums.resize(w * h);
		flux.resize(w * h);
		hitPoints.resize(w * h);
		directillum.resize(w * h);
		c.resize(w * h);
		for (int i = 0; i < w * h; ++i) {
			radius2[i] = initialRadius * initialRadius;
			photonNums[i] = 0;
			flux[i] = Vec();
		}
	}

	HashGrid hashGrid;
	std::shared_ptr<Sampler> sampler;
	std::vector<double> radius2;
	std::vector<long long> photonNums;
	std::vector<Vec> flux;
	std::vector<Vec> directillum;
	std::vector<HPoint> hitPoints;
	std::vector<Vec> c;
	const double initialRadius;
	const int maxDepth;
	const int nIterations;
	const long long nPhotonsPerRenderStage;
	const double alpha;
};


class Renderer {
public:
	Renderer(const std::shared_ptr<Scene> &pScene, const std::shared_ptr<Integrator> &pIntegrator,
		const std::shared_ptr<Film> &pFilm): scene(pScene), integrator(pIntegrator), film(pFilm) {

	}
	void Render() {
		integrator->Render(*scene);
		film->SaveImage();
	}
private:
	std::shared_ptr<Scene> scene;
	std::shared_ptr<Integrator> integrator;
	std::shared_ptr<Film> film;
};

int main(int argc, char *argv[]) {
	
	clock_t begin = clock();

	int w = 1024, h = 768;
	int nIterations = (argc == 2) ? atol(argv[1]) : 256; //(argc == 2) ? std::max(atoll(argv[1]) / render_stage_number, (long long)1) : render_stage_number;

	//std::cout << nIterations << std::endl;

	//Initialize(w, h, 0.8);
	//hashGridSpinlocks.resize(w * h);

	std::shared_ptr<Film> film = std::shared_ptr<Film>(new Film(w, h));

	std::shared_ptr<Scene> scene = std::shared_ptr<Scene>(new Scene);

	// trace eye rays and store measurement points
	//Ray cam(Vec(50, 48, 295.6), Vec(0, -0.042612, -1).norm());
	//Vec cx = Vec(w*.5135 / h), cy = (cx%cam.d).norm()*.5135, *c = new Vec[w*h], vw;

	Vec camPos(50, 52, 295.6), cz(0, -0.042612, -1);
	//double filmDis = cz.length();
	double filmDis = 1.0;
	Vec cx = Vec(w * .5135 / h).norm();
	Vec cy = (cx % cz).norm();
	double fovy = 28.7993;

	//std::cout << camPos + cz << std::endl << cy << std::endl;
	
	std::shared_ptr<Camera> camera = std::shared_ptr<Camera>(new PinHoleCamera(film, camPos, cz, cx, cy, fovy, filmDis));
	std::shared_ptr<Sampler> randomSampler = std::shared_ptr<Sampler>(new RandomSampler(123));
	std::shared_ptr<Sampler> haltonSampler = std::shared_ptr<Sampler>(new HaltonSampler());
	std::shared_ptr<Integrator> integrator = std::shared_ptr<Integrator>(new SPPM(nIterations, render_stage_number, 20, 0.8, haltonSampler));
	fprintf(stderr, "Load Scene ...\n");
	//scene->SetCamera(cam, cx, cy);
	scene->SetCamera(camera);
	scene->AddShape(std::shared_ptr<Shape>(new Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(), Vec(.75, .25, .25), DIFF)));//Left
	scene->AddShape(std::shared_ptr<Shape>(new Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(), Vec(.25, .25, .75), DIFF)));//Right
	scene->AddShape(std::shared_ptr<Shape>(new Sphere(1e5, Vec(50, 40.8, 1e5), Vec(), Vec(.75, .75, .75), DIFF)));//Back
	//scene->AddShape(std::shared_ptr<Shape>(new Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Vec(), Vec(.75, .75, .75), DIFF)));//Frnt
	scene->AddShape(std::shared_ptr<Shape>(new Sphere(1e5, Vec(50, 1e5, 81.6), Vec(), Vec(.75, .75, .75), DIFF)));//Botm
	scene->AddShape(std::shared_ptr<Shape>(new Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(), Vec(.75, .75, .75), DIFF)));//Top
	scene->AddShape(std::shared_ptr<Shape>(new Sphere(16.5, Vec(27, 16.5, 47), Vec(), Vec(1, 1, 1)*.999, REFR)));//Mirr
	scene->AddShape(std::shared_ptr<Shape>(new Sphere(7.0, Vec(27, 16.5, 47), Vec(), Vec(.25, .25, .75), DIFF)));//Mirr
	scene->AddShape(std::shared_ptr<Shape>(new Sphere(16.5, Vec(73, 26.5, 78), Vec(), Vec(1, 1, 1)*.999, REFR)));//Glass
	scene->AddShape(std::shared_ptr<Shape>(new Sphere(9.5, Vec(53, 9.5, 88), Vec(), Vec(1, 1, 1)*.999, REFR)));//Glass
	std::shared_ptr<Shape> lightShape = std::shared_ptr<Shape>(new Sphere(8.0, Vec(50, 81.6 - 16.5, 81.6), Vec(0.3, 0.3, 0.3) * 100, Vec(), DIFF));//Lite
	std::shared_ptr<Light> light0 = std::shared_ptr<Light>(new AreaLight(lightShape));
	scene->AddLight(light0);
	scene->Initialize();
	film->SetFileName("cornellbox14.bmp");
	std::shared_ptr<Renderer> renderer = std::shared_ptr<Renderer>(new Renderer(scene, integrator, film));
	renderer->Render();

	clock_t end = clock();

	std::cout << "cost time: "<< (end - begin) / 1000.0 / 60.0 <<" min"<< std::endl;



	//test
	/*
	Intersection isect;
	isect.n = Vec(0, 1, 0);
	isect.nl = isect.n * -1;
	isect.wo = Vec(1, -1, 0) * -1;
	Vec wi;
	double pdf;
	TransmissionBSDF bsdf(isect);
	bsdf.Sample_f(isect.wo, &wi, &pdf, Vec(0.5, 0.5, 0.5));
	std::cout << "------------------------" << std::endl;

	bool into = isect.n.dot(isect.nl) > 0;                // Ray from outside going in?
	double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = (isect.wo * -1).dot(isect.nl), cos2t;
	cos2t = 1 - nnt * nnt*(1 - ddn * ddn);
	Vec tdir = ((isect.wo * -1)*nnt - isect.n * ((into ? 1 : -1)*(ddn*nnt + sqrt(cos2t)))).norm();
	double a = nt - nc, b = nt + nc, R0 = a * a / (b*b), c = 1 - (into ? -ddn : tdir.dot(isect.n));
	double Re = R0 + (1 - R0)*c*c*c*c*c;
	std::cout << tdir << std::endl << std::endl << Re << std::endl;*/
	/**
	// save the image after tone mapping and gamma correction
	FILE* f = fopen("cornellbox10.png", "wb");
	unsigned char *RGBs = new unsigned char[w * h * 3];
	for (int j = 0; j < h; ++j) {
		for (int i = 0; i < w; ++i) {
			int index = 3 * j * w + 3 * i;
			int index_c = j * w + i;
			RGBs[index] = (unsigned char)(toInt(c[index_c].x));
			RGBs[index + 1] = (unsigned char)(toInt(c[index_c].y));
			RGBs[index + 2] = (unsigned char)(toInt(c[index_c].z));
		}
	}
	svpng(f, w, h, RGBs, 0);
	delete[] RGBs;
	fclose(f);*/
	//_CrtDumpMemoryLeaks();
}
