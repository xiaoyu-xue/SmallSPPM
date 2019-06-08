#include "svpng.inc"
#include <math.h>  
#include <stdlib.h> 
#include <stdio.h>  
#include <Random>
#include <vector>
#include <crtdbg.h>
#include <iostream>
#include <algorithm>
#include <atomic>
#include <ctime>
#define _CRTDBG_MAP_ALLOC

const double PI = 3.14159265358979;
const double INV_PI = 0.31830988618379067154;
const double ALPHA = 0.7;
const int render_stage_number = 7000000;
const double PiOver2 = 1.57079632679489661923;
const double PiOver4 = 0.78539816339744830961;
const double eps = 10e-6;
const double Inf = 1e20;


class Shape;
class BSDF;


std::mt19937_64 rng(1234);
std::uniform_real_distribution<double> uniform;
enum ReflectionType { DIFF, SPEC, REFR };  // material types, used in radiance()

double Random() {
	return uniform(rng);
}

// Halton sequence with reverse permutation
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
}



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
	inline Vec norm() { return (*this) * (1.0 / sqrt(x*x + y * y + z * z)); }
	inline double dot(const Vec &b) const { return x * b.x + y * b.y + z * b.z; }
	Vec operator%(Vec&b) const { return Vec(y*b.z - z * b.y, z*b.x - x * b.z, x*b.y - y * b.x); }
	double& operator[](int i) { return i == 0 ? x : i == 1 ? y : z; }
	double maxValue() const {
		return std::max(x, std::max(y, z));
	}
};

Vec operator*(double a, Vec b) { return Vec(a * b.x, a * b.y, a * b.z); }

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
			return Fa / cosTheta;
		}
		Vec td = ((-1 * wo) * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
		double Re = Fresnell(wo, td, n, nl);
		double P = Re * 0.5 + 0.25;
		if (rand.z < P) {
			*wi = (nl * 2.0 * nl.dot(wo) - wo).norm();
			*pdf = 1.0;
			double cosTheta = std::abs((*wi).dot(n));
			return Fa * Re / cosTheta / P;
		}
		else {

			*wi = td;
			*pdf = 1.0;
			double cosTheta = std::abs((*wi).dot(n));
			return Fa * (1.0 - Re) / cosTheta / (1 - P);
		}
	}

	Vec f(const Vec &wo, const Vec &wi) const {
		return Vec();
	}

	double Fresnell(const Vec &wo, const Vec &td, const Vec &n, const Vec &nl) const {
		bool into = (n.dot(nl) > 0.0);
		double nnt = into ? nc / nt : nt / nc, ddn = (-1 * wo).dot(nl), cos2t;
		cos2t = 1 - nnt * nnt*(1 - ddn * ddn);
		double a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : td.dot(n));
		double Re = R0 + (1 - R0) * c * c* c * c * c;
		return Re;
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
	virtual Vec SampleLight(Vec *pos, Vec *dir, Vec *lightNorm, double *pdfPos, double *pdfDir) const {
		SampleOnLight(pos, dir, lightNorm, pdfPos, pdfDir);
		return Emission();
	}
	virtual int GetId() const = 0;
	virtual bool IsAreaLight() const { return false; }
	virtual std::shared_ptr<Shape> GetShapePtr() const = 0;
protected:
	virtual void SampleOnLight(Vec *pos, Vec *dir, Vec *lightNorm, double *pdfPos, double *pdfDir) const = 0;
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

	bool IsAreaLight() const { return true; }

	std::shared_ptr<Shape> GetShapePtr() const { return shape; }
protected:
	void SampleOnLight(Vec *pos, Vec *dir, Vec *lightNorm, double *pdfPos, double *pdfDir) const {
		//sample a position
		*pos = shape->Sample(pdfPos, Vec(Random(), Random(), Random()));
		*lightNorm = shape->GetNorm(*pos);
		Vec ss, ts;
		CoordinateSystem(*lightNorm, &ss, &ts);
		Vec dirLocal = CosineSampleHemisphere(Vec(Random(), Random(), Random()));
		double cosTheta = dirLocal.z;
		*dir = (ss * dirLocal.x + ts * dirLocal.y + *lightNorm * dirLocal.z).norm();
		*pdfDir = CosineHemispherePdf(cosTheta);
	}
private:
	std::shared_ptr<Shape> shape;
};



std::vector<double> radius2;
std::vector<long long> photonNums;
std::vector<Vec> flux;
std::vector<Vec> directillum;
std::vector<HPoint> hitPoints;
std::vector<std::vector<HPoint*>> hashGrid;
long long hashNum, pixelIndex;
double hashCellSize;
AABB hpbbox;





void ClearHashGrid() {
	for (auto &e : hashGrid) {
		std::vector<HPoint*>().swap(e);
	}
}

// spatial hash function
inline unsigned int hash(const int ix, const int iy, const int iz) {
	return (unsigned int)((ix * 73856093) ^ (iy * 19349663) ^ (iz * 83492791)) % hashNum;
}


void BuildHashGrid(const int w, const int h, int stage_num) {
	// find the bounding box of all the measurement points
	hpbbox.reset();
	for (HPoint &hp :hitPoints) {
		if (hp.used) hpbbox.fit(hp.pos);
	}

	// heuristic for initial radius
	Vec ssize = hpbbox.maxPoint - hpbbox.minPoint;
	double irad = ((ssize.x + ssize.y + ssize.z) / 3.0) / ((w + h) / 2.0) * 2.0;
	// determine hash table size
	// we now find the bounding box of all the measurement points inflated by the initial radius
	hpbbox.reset();
	int vphoton = 0;

	for (HPoint &hp : hitPoints) {
		if (!hp.used) continue;
		if (stage_num == 0) {
			radius2[hp.pix] = irad * irad;
			photonNums[hp.pix] = 0;
			flux[hp.pix] = Vec();
		}
		vphoton++;
		hpbbox.fit(hp.pos - irad);
		hpbbox.fit(hp.pos + irad);
	}

	// make each grid cell two times larger than the initial radius
	hashCellSize = 1.0 / (irad*2.0);
	hashNum = vphoton;

	// build the hash table

	hashGrid.resize(hashNum);


	for (HPoint &hp : hitPoints) {
		if (!hp.used) continue;
		Vec BMin = ((hp.pos - irad) - hpbbox.minPoint) * hashCellSize;
		Vec BMax = ((hp.pos + irad) - hpbbox.minPoint) * hashCellSize;
		for (int iz = abs(int(BMin.z)); iz <= abs(int(BMax.z)); iz++)
		{
			for (int iy = abs(int(BMin.y)); iy <= abs(int(BMax.y)); iy++)
			{
				for (int ix = abs(int(BMin.x)); ix <= abs(int(BMax.x)); ix++)
				{
					int hv = hash(ix, iy, iz);
					hashGrid[hv].push_back(&hp);
				}
			}
		}
	}
}


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

	void SetCamera(const Ray &camera, const Vec &cX, const Vec &cY) {
		cam = camera;
		cx = cX;
		cy = cY;
		shapeNum = 0;
	}

	void AddShape(std::shared_ptr<Shape> shape) {
		shapes.push_back(shape);
		shapes[shapeNum]->shapeId = shapeNum;
		++shapeNum;
	}

	void AddLight(std::shared_ptr<Light> light) {
		lights.push_back(light);
		AddShape(light->GetShapePtr());
	}

	bool Intersect(const Ray &r, double *t, Intersection *isect, std::shared_ptr<Shape> &hitObj) const {
		int n = shapes.size();
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
private:
	std::vector<std::shared_ptr<Shape>> shapes;
	std::vector<std::shared_ptr<Light>> lights;
	Ray cam;
	Vec cx, cy;
	int shapeNum;
};

// tone mapping and gamma correction
int toInt(double x) {
	return int(pow(1 - exp(-x), 1 / 2.2) * 255 + .5);
}


void GenratePhoton(const Scene &scene, Ray *pr, Vec *f) {
	const std::vector<std::shared_ptr<Light>> lights = scene.GetLights();
	int lightsNum = lights.size();
	double lightPdf = 1.0 / lightsNum;
	int lightindex = (int)(Random() * lightsNum);
	const Light &light = *(lights[lightindex]);
	Vec Le = light.Emission();
	Vec pos, lightDir, lightNorm;
	double pdfPos, pdfDir;
	light.SampleLight(&pos, &lightDir, &lightNorm, &pdfPos, &pdfDir);
	pr->o = pos + lightDir * eps;
	pr->d = lightDir;
	double cosTheta = std::abs(lightNorm.dot(lightDir));
	*f = Le * cosTheta / (pdfPos * pdfDir);
}



Vec DirectIllumination(const Scene &scene, const Intersection &isect, const std::shared_ptr<BSDF> &bsdf, Vec importance, Vec u) {
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


void TraceEyePath(const Scene &scene, const Ray &ray, int maxDepth, long long pixel) {
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
			Vec f = bsdf->Sample_f(-1 * r.d, &wi, &pdf, Vec(Random(), Random(), Random()));
			importance = f * std::abs(wi.dot(isect.n)) * importance / pdf;
			r.o = isect.hit;
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
			hitPoints[pixel] = hp;
			if ((i == 0 || deltaBoundEvent) && hitObj->IsLight())
				directillum[hp.pix] = directillum[hp.pix] + importance * hitObj->GetEmission();
			else
				directillum[hp.pix] = directillum[hp.pix] + 
					DirectIllumination(scene, isect, bsdf, hp.importance, Vec(Random(), Random(), Random()));

			return;
		}
	}
}

void TracePhoton(const Scene &scene, const Ray &ray, Vec photonFlux, int maxDepth) {
	Ray r = ray;
	for (int i = 0; i < maxDepth; ++i) {
		double t;
		Intersection isect;
		std::shared_ptr<Shape> hitObj;
		if (!scene.Intersect(r, &t, &isect, hitObj)) return;
		std::shared_ptr<BSDF> bsdf = hitObj->GetBSDF(isect);
		Vec wi;
		double pdf;
		Vec f = bsdf->Sample_f(-1 * r.d, &wi, &pdf, Vec(Random(), Random(), Random()));
		Vec estimation = f * std::abs(wi.dot(isect.n)) / pdf;
		if (bsdf->IsDelta()) {
			photonFlux = photonFlux * estimation;
			r.o = isect.hit;
			r.d = wi;
		}
		else {
			if (i > 0) {
				// photon ray
				// find neighboring measurement points and accumulate flux via progressive density estimation
				Vec hh = (isect.hit - hpbbox.minPoint) * hashCellSize;
				int ix = std::abs(int(hh.x)), iy = std::abs(int(hh.y)), iz = std::abs(int(hh.z));
				// strictly speaking, we should use #pragma omp critical here.
				// it usually works without an artifact due to the fact that photons are 
				// rarely accumulated to the same measurement points at the same time (especially with QMC).
				// it is also significantly faster.
				std::vector<HPoint*> &hp = hashGrid[hash(ix, iy, iz)];
				for (HPoint* hitpoint : hp) {
					Vec v = hitpoint->pos - isect.hit;
					if ((hitpoint->nrm.dot(isect.n) > 1e-3) && (v.dot(v) <= radius2[hitpoint->pix])) {
						// unlike N in the paper, hitpoint->n stores "N / ALPHA" to make it an integer value
						double g = (photonNums[hitpoint->pix] * ALPHA + ALPHA) / (photonNums[hitpoint->pix] * ALPHA + 1.0);
						radius2[hitpoint->pix] = radius2[hitpoint->pix] * g;
						photonNums[hitpoint->pix]++;
						Vec contribution = hitpoint->importance * bsdf->f(hitpoint->outDir, -1 * r.d) * photonFlux;
						flux[hitpoint->pix] = (flux[hitpoint->pix] + contribution) * g;
					}
				}
			}

			double p = estimation.maxValue();
			if (p < 1) {
				if (Random() < p) photonFlux = photonFlux / p;
				else break;
			}
			photonFlux = photonFlux * estimation;
			r.o = isect.hit;
			r.d = wi;
		}
	}
}



int main(int argc, char *argv[]) {
	
	clock_t begin = clock();

	int w = 1024 * 2, h = 768 * 2;
	long long nIterations = (argc == 2) ? atoll(argv[1]) : 256; //(argc == 2) ? std::max(atoll(argv[1]) / render_stage_number, (long long)1) : render_stage_number;

	std::cout << nIterations << std::endl;
	radius2.resize(w * h);
	photonNums.resize(w * h);
	flux.resize(w * h);
	hitPoints.resize(w * h);
	directillum.resize(w * h);

	std::shared_ptr<Scene> scene = std::shared_ptr<Scene>(new Scene);

	// trace eye rays and store measurement points
	Ray cam(Vec(50, 48, 295.6), Vec(0, -0.042612, -1).norm());
	Vec cx = Vec(w*.5135 / h), cy = (cx%cam.d).norm()*.5135, *c = new Vec[w*h], vw;

	fprintf(stderr, "Load Scene ...\n");
	scene->SetCamera(cam, cx, cy);
	scene->AddShape(std::shared_ptr<Shape>(new Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(), Vec(.75, .25, .25), DIFF)));
	scene->AddShape(std::shared_ptr<Shape>(new Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(), Vec(.25, .25, .75), DIFF)));
	scene->AddShape(std::shared_ptr<Shape>(new Sphere(1e5, Vec(50, 40.8, 1e5), Vec(), Vec(.75, .75, .75), DIFF)));//Back
	scene->AddShape(std::shared_ptr<Shape>(new Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Vec(), Vec(.75, .75, .75), DIFF)));//Frnt
	scene->AddShape(std::shared_ptr<Shape>(new Sphere(1e5, Vec(50, 1e5, 81.6), Vec(), Vec(.75, .75, .75), DIFF)));//Botm
	scene->AddShape(std::shared_ptr<Shape>(new Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(), Vec(.75, .75, .75), DIFF)));//Top
	scene->AddShape(std::shared_ptr<Shape>(new Sphere(16.5, Vec(27, 16.5, 47), Vec(), Vec(1, 1, 1)*.999, REFR)));//Mirr
	scene->AddShape(std::shared_ptr<Shape>(new Sphere(7.0, Vec(27, 16.5, 47), Vec(), Vec(.25, .25, .75), DIFF)));//Mirr
	scene->AddShape(std::shared_ptr<Shape>(new Sphere(16.5, Vec(73, 26.5, 78), Vec(), Vec(1, 1, 1)*.999, REFR)));//Glas
	scene->AddShape(std::shared_ptr<Shape>(new Sphere(9.5, Vec(53, 9.5, 88), Vec(), Vec(1, 1, 1)*.999, REFR)));
	std::shared_ptr<Shape> lightShape = std::shared_ptr<Shape>(new Sphere(8.0, Vec(50, 81.6 - 16.5, 81.6), Vec(0.3, 0.3, 0.3) * 100, Vec(), DIFF));
	std::shared_ptr<Light> light0 = std::shared_ptr<Light>(new AreaLight(lightShape));
	scene->AddLight(light0);

	fprintf(stderr, "Rendering ...\n");
	for (int render_stage = 0; render_stage < nIterations; ++render_stage) {
#pragma omp parallel for schedule(guided)
		for (int y = 0; y < h; y++) {
			//fprintf(stderr, "\rHitPointPass %5.2f%%", 100.0*y / (h - 1));
//#pragma omp parallel for schedule(dynamic, 1)
			for (int x = 0; x < w; x++) {
				int pixel = x + y * w;
				hitPoints[pixel].used = false;
				double u = Random() - 0.5, v = Random() - 0.5;
				Vec d = cx * ((x + 0.5 + u) / w - 0.5) + cy * (-(y + 0.5 + v) / h + 0.5) + cam.d;
				//trace(Ray(cam.o + d * 140, d.norm()), 0, true, Vec(), Vec(1, 1, 1), 0, false, pixel);
				Ray ray(cam.o + d * 140, d.norm());
				TraceEyePath(*scene, ray, 21, pixel);
			}
		}

		BuildHashGrid(w, h, render_stage);
		{
			//fprintf(stderr, "\n");
			double percentage = 100.*(render_stage + 1) / nIterations;
			//fprintf(stderr, "\rPhotonPass %5.2f%%", percentage);
			Ray ray;
			Vec photonFlux;
#pragma omp parallel for schedule(guided)
			for (int j = 0; j < render_stage_number; j++) {
				GenratePhoton(*scene, &ray, &photonFlux);
				TracePhoton(*scene, ray, photonFlux, 21);
				//GenratePhoton(&r, &f, m + j);
				//trace(r, 0, 0 > 1, f, vw, m + j);
			}
			//fprintf(stderr, "\n");
			fprintf(stderr, "\rPhotonPass %5.2f%%", percentage);
		}
		ClearHashGrid();
	}

	// density estimation
	
	//for (int i = 0; i < flux.size(); ++i) {
	//	std::cout << flux[i].x << " " << flux[i].y << " " << flux[i].z << std::endl;
	//}
	std::cout << "\nflux size: " << flux.size() << std::endl;
#pragma omp parallel for schedule(guided)
	for (int i = 0; i < flux.size(); ++i) {
		c[i] = c[i] + flux[i] * (1.0 / (PI * radius2[i] * nIterations * render_stage_number))
			+ directillum[i] / nIterations;
		//c[i] = c[i] + directillum[i] / nIterations;
		//std::cout << c[i].x << " " << c[i].y << " " << c[i].z << std::endl;
	}

	clock_t end = clock();
	std::cout << "cost time "<< (end - begin) / 1000.0 / 60.0 <<"min"<< std::endl;
	// save the image after tone mapping and gamma correction
	FILE* f = fopen("cornellbox9.png", "wb");
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
	fclose(f);
	//_CrtDumpMemoryLeaks();
}
