#include "svpng.inc"
#include <math.h>  
#include <stdlib.h> 
#include <stdio.h>  
#include <Random>
#include <vector>
#include <crtdbg.h>
#include <iostream>
#include <algorithm>
#define _CRTDBG_MAP_ALLOC

const double PI = 3.14159265358979;
const double INV_PI = 0.31830988618379067154;
const double ALPHA = 0.7;
const int render_stage_number = 200000;
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
	BSDF(const Vec &normal, const Vec &normalL) : n(normal), nl(normalL) {}
	virtual double Pdf(const Vec &wo, const Vec &wi) const = 0;
	virtual Vec Sample_f(const Vec &wo, Vec *wi, double *pdf, Vec rand = Vec(0, 0, 0)) const = 0;
	virtual Vec f(const Vec &wo, const Vec &wi) const { return Vec(0, 0, 0); }
	virtual bool IsDelta() const { return false; }
protected:
	const Vec n, nl;
};

class DiffuseBSDF : public BSDF {
public:
	DiffuseBSDF(const Vec &n, const Vec &nl, Vec r) : BSDF(n, nl), R(r) {}

	double Pdf(const Vec &wo, const Vec &wi) const {
		return std::abs(wi.dot(nl)) * INV_PI;
	}

	Vec Sample_f(const Vec &wo, Vec *wi, double *pdf, Vec rand) const {
		double r1 = 2. * PI * rand[0], r2 = rand[1];
		double r2s = sqrt(r2);
		Vec w = nl, u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).norm();
		Vec v = w % u;
		*wi = (u*cos(r1)*r2s + v * sin(r1)*r2s + w * sqrt(1 - r2)).norm();
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
	SpecularBSDF(const Vec &n, const Vec &nl, Vec r = Vec(1.0, 1.0, 1.0)) : BSDF(n, nl), R(r) {}

	double Pdf(const Vec &wo, const Vec &wi) const {
		return 1.0;
	}

	Vec Sample_f(const Vec &wo, Vec *wi, double *pdf, Vec rand) const {
		*wi = nl * 2.0 * nl.dot(wo) - wo;
		*pdf = Pdf(wo, *wi);
		return f(wo, *wi);
	}

	Vec f(const Vec &wo, const Vec &wi) const {
		return R;
	}

	bool IsDelta() const { return true; }
private:
	Vec R;
};

class TransmissionBSDF : public BSDF {
public:
	TransmissionBSDF(const Vec &n, const Vec &nl, Vec fa = Vec(1.0, 1.0, 1.0), double eta1 = 1.0, double eta2 = 1.5) :
		BSDF(n, nl), Fa(fa), nc(eta1), nt(eta2) {}

	double Pdf(const Vec &wo, const Vec &wi) const {
		return 1.0;
	}

	Vec Sample_f(const Vec &wo, Vec *wi, double *pdf, Vec rand) const {
		bool into = (n.dot(nl) > 0.0);
		double nnt = into ? nc / nt : nt / nc, ddn = (-1 * wo).dot(nl), cos2t;
		// total internal reflection
		if ((cos2t = 1 - nnt * nnt*(1 - ddn * ddn)) < 0) {
			*wi = nl * 2.0 * nl.dot(wo) - wo;
			*pdf = Pdf(wo, *wi);
			return Fa;
		}
		Vec td = ((-1 * wo) * nnt - nl * (ddn * nnt + sqrt(cos2t))).norm();
		double Re = Fresnell(wo, nl, into);
		double P = Re * 0.5 + 0.25;
		if (P < rand[2]) {
			*wi = nl * 2.0 * nl.dot(wo) - wo;
			*pdf = Pdf(wo, *wi);
			return Fa * Re / P;
		}
		else {
			*wi = td;
			*pdf = Pdf(wo, *wi);
			return Fa * (1.0 - Re) / (1 - P);
		}
	}

	Vec f(const Vec &wo, const Vec &wi) const {
		return Fa;
	}

	bool IsDelta() const { return true; }

private:
	double Fresnell(const Vec &wo, const Vec &nl, bool into) const {
		double nnt = into ? nc / nt : nt / nc, ddn = (-1 * wo).dot(nl), cos2t;
		cos2t = 1 - nnt * nnt*(1 - ddn * ddn);
		Vec td = ((-1 * wo) * nnt - nl * (ddn * nnt + sqrt(cos2t))).norm();
		double a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : std::abs(td.dot(nl)));
		double Re = R0 + (1 - R0) * c * c* c * c * c;
		return Re;
	}

	double nc, nt;
	Vec Fa;
};

class Shape {
public:
	Shape(ReflectionType type, const Vec &color, bool isL = false): reflType(type), c(color), isLight(isL){}
	virtual double Intersect(const Ray &r) const = 0;
	virtual Vec Sample(Vec rand) const = 0;
	virtual std::shared_ptr<BSDF> GetBSDF(const Vec &n, const Vec &nl) const {
		if (reflType == DIFF) {
			return std::shared_ptr<BSDF>(new DiffuseBSDF(n, nl, c));
		}
		else if (reflType == SPEC) {
			return std::shared_ptr<BSDF>(new SpecularBSDF(n, nl, c));
		}
		else if (reflType == REFR) {
			return std::shared_ptr<BSDF>(new TransmissionBSDF(n, nl, c));
		}
	}
	bool IsLight() const { return isLight; }
	virtual Vec GetNorm(const Vec &point) const = 0;
	virtual Vec GetEmission() const = 0;
	int GetId() const { return shapeId; }
	friend class Scene;
private:
	ReflectionType reflType;
	Vec c;
	bool isLight;
	int shapeId;

};

class Light {
public:
	Light() {}
	virtual Vec DirectIllumination(const Vec &hitpoint, const Vec &wo, const std::shared_ptr<BSDF> bsdf, 
		const Vec &nl, const Vec &importance, Vec *dir, Vec u) const = 0;
	virtual Vec Emission() const = 0;
	virtual Vec SampleLight(Vec *pos, Vec *dir, Vec *lightNorm, double *pdfPos, double *pdfDir) const {
		SampleOnLight(pos, dir, lightNorm, pdfPos, pdfDir);
		return Emission();
	}
	virtual int GetId() const = 0;
	friend class Scene;
protected:
	virtual void SampleOnLight(Vec *pos, Vec *dir, Vec *lightNorm, double *pdfPos, double *pdfDir) const = 0;
};


class Sphere: public Shape {
public:
	Sphere(double r_, Vec p_, Vec e_, Vec c_, ReflectionType reflType) : 
		rad(r_), p(p_), e(e_), c(c_), Shape(reflType, c_){}

	double Intersect(const Ray &r) const {
		// ray-sphere Intersection returns distance
		Vec op = p - r.o;
		double t, b = op.dot(r.d), det = b * b - op.dot(op) + rad * rad;
		if (det < 0) {
			return 1e20;
		}
		else {
			det = sqrt(det);
		}
		return (t = b - det) > 1e-4 ? t : ((t = b + det) > 1e-4 ? t : 1e20);
	}

	Vec Sample(Vec rand) const {
		return UniformSampleSphere(rand);
	}

	Vec Emission() const {
		return e;
	}

	Vec GetNorm(const Vec & point) const {
		return (point - p).norm();
	}

	Vec GetEmission() const {
		return e;
	}
double rad; Vec p, e, c;

protected:
	
};


class SphereLight : public Light{
public:
	SphereLight(const std::shared_ptr<Sphere>& sph): sphere(sph){}
	Vec DirectIllumination(const Vec &hitpoint, const Vec &wo, const std::shared_ptr<BSDF> bsdf,
		const Vec &nl, const Vec &importance, Vec *dir, Vec u) const {
		Vec sw = sphere->p - hitpoint, su = ((fabs(sw.x) > .1 ? Vec(0, 1) : Vec(1)) % sw).norm(), sv = sw % su;
		double cos_a_max = sqrt(1 - sphere->rad * sphere->rad / (hitpoint - sphere->p).dot(hitpoint - sphere->p));
		double zeta1 = u.x, zeta2 = u.y;
		double cos_a = 1 - zeta1 + zeta1 * cos_a_max;
		double sin_a = sqrt(1 - cos_a * cos_a);
		double phi = 2 * PI * zeta2;
		Vec l = su * cos(phi) * sin_a + sv * sin(phi) * sin_a + sw * cos_a;
		*dir = l.norm();
		double omega = 2 * PI*(1 - cos_a_max);
		Vec f = bsdf->f(wo, *dir);
		return importance * f * l.dot(nl) * sphere->e * omega;  // 1/pi for brdf
	}

	Vec Emission() const {
		return sphere->e;
	}

	int GetId() const {
		return sphere->GetId();
	}
protected:
	void SampleOnLight(Vec *pos, Vec *dir, Vec *lightNorm, double *pdfPos, double *pdfDir) const {
		//sample a position
		*pos = UniformSampleSphere(Vec(Random(), Random(), Random())) * sphere->rad + sphere->p;
		*pdfPos = 1.f / (4.0 * PI * sphere->rad * sphere->rad);
		*lightNorm = (*pos - sphere->p).norm();
		Vec ss, ts;
		CoordinateSystem(*lightNorm, &ss, &ts);
		Vec dirLocal = CosineSampleHemisphere(Vec(Random(), Random(), Random()));
		double cosTheta = dirLocal.z;
		Vec lightDir = (ss * dirLocal.x + ts * dirLocal.y + *lightNorm * dirLocal.z).norm();
		*pdfDir = CosineHemispherePdf(cosTheta);
	}
private:
	std::shared_ptr<Sphere> sphere;
};

std::vector<double> radius2;
std::vector<unsigned int> photonNums;
std::vector<Vec> flux;
std::vector<Vec> directillum;
std::vector<HPoint> hitPoints;
std::vector<std::vector<HPoint*>> hashGrid;
unsigned int hashNum, pixelIndex, photonTotalNum;
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

	void AddLight(std::shared_ptr<Light> light, std::shared_ptr<Shape> shape) {
		lights.push_back(light);
		AddShape(shape);
	}

	bool Intersect(const Ray &r, double &t, std::shared_ptr<Shape> &hitObj) const {
		int n = shapes.size();
		double d;
		t = Inf;
		for (int i = 0; i < n; ++i) {
			d = shapes[i]->Intersect(r);
			if (d < t) {
				t = d;
				hitObj = shapes[i];
			}
		}
		return t < Inf;
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

// find the closet interection
inline bool Intersect(const Ray &r, double &t, int &id) {
	int n = sizeof(sph) / sizeof(Sphere);
	double d, inf = 1e20; t = inf;
	for (int i = 0; i < n; i++) {
		d = sph[i].Intersect(r);
		if (d < t) {
			t = d;
			id = i;
		}
	}
	return t < inf;
}



void GenratePhoton(const Scene &scene, Ray *pr, Vec *f) {
	const std::vector<std::shared_ptr<Light>> lights = scene.GetLights();
	int lightsNum = lights.size();
	double lightPdf = 1.0 / lightsNum;
	int lightindex = (int)(Random() * lightsNum);
	const Light &light = *lights[lightindex];
	Vec Le = light.Emission();
	Vec pos, lightDir, lightNorm;
	double pdfPos, pdfDir;
	light.SampleLight(&pos, &lightDir, &lightNorm, &pdfPos, &pdfDir);
	pr->o = pos + lightDir * eps;
	pr->d = lightDir;
	double cosTheta = std::abs(lightNorm.dot(lightDir));
	*f = Le * cosTheta / (pdfPos * pdfDir);
}

/*
void GenratePhoton(Ray *pr, Vec *f, int i) {
	int lightindex = (int)(sizeof(sph) / sizeof(Sphere)) - 1;
	const Sphere &light = sph[lightindex];
	Vec Le = light.e;
	//sample a position
	Vec pos = UniformSampleSphere(Vec(Random(), Random(), Random())) * light.rad + light.p;
	double Pdfpos = 1.f / (4 * PI * light.rad * light.rad);
	Vec lightNorm = (pos - light.p).norm();
	Vec ss, ts;
	CoordinateSystem(lightNorm, &ss, &ts);
	Vec dirLocal = CosineSampleHemisphere(Vec(Random(), Random(), Random()));
	double cosTheta = dirLocal.z;
	Vec lightDir = (ss * dirLocal.x + ts * dirLocal.y + lightNorm * dirLocal.z).norm();
	double Pdfdir = CosineHemispherePdf(cosTheta);
	pr->o = pos + lightDir * eps;
	pr->d = lightDir;
	*f = Le * cosTheta / (Pdfpos * Pdfdir);
}*/


Vec DirectIllumination(const Scene &scene, Vec hitpoint, Vec wo, const std::shared_ptr<BSDF> bsdf, Vec nl, Vec importance, Vec u) {
	Vec L;
	const std::vector<std::shared_ptr<Light>> lights = scene.GetLights();
	for (auto light : lights) {
		Vec dir;
		std::shared_ptr<Shape> hitObj;
		double t;
		Vec Li = light->DirectIllumination(hitpoint, wo, bsdf, nl, importance, &dir, u);
		if (scene.Intersect(Ray(hitpoint, dir), t, hitObj) && hitObj == light) {
			L = L + Li;
		}
	}
	return L;
}


Vec DirectIllumination(Vec hitpoint, Vec nl, Vec f, Vec u) {
	Vec e;
	int numSpheres = (int)(sizeof(sph) / sizeof(Sphere));
	for (int i = 0; i < numSpheres; i++) {
		const Sphere &s = sph[i];
		if (s.e.x <= 0 && s.e.y <= 0 && s.e.z <= 0) continue; // skip non-lights

		Vec sw = s.p - hitpoint, su = ((fabs(sw.x) > .1 ? Vec(0, 1) : Vec(1)) % sw).norm(), sv = sw % su;
		double cos_a_max = sqrt(1 - s.rad*s.rad / (hitpoint - s.p).dot(hitpoint - s.p));
		double eps1 = u.x, eps2 = u.y;
		double cos_a = 1 - eps1 + eps1 * cos_a_max;
		double sin_a = sqrt(1 - cos_a * cos_a);
		double phi = 2 * PI*eps2;
		Vec l = su * cos(phi)*sin_a + sv * sin(phi)*sin_a + sw * cos_a;
		l = l.norm();
		double t;
		int id;
		if (Intersect(Ray(hitpoint, l), t, id) && id == i) {  // shadow ray
			double omega = 2 * PI*(1 - cos_a_max);
			e = e + f.mul(s.e*l.dot(nl)*omega) * INV_PI;  // 1/pi for brdf
		}
	}
	return e;
}

/*
void trace(const Ray &r, int dpt, bool m, const Vec &fl, const Vec &adj, int i, bool specularbounce = false, int pixel = -1)
{
	double t;
	int id;

	dpt++;
	if (!Intersect(r, t, id) || (dpt >= 20))return;

	int d3 = dpt * 3;
	const Sphere &obj = sph[id];
	Vec x = r.o + r.d*t, n = (x - obj.p).norm(), f = obj.c;
	Vec nl = n.dot(r.d) < 0 ? n : n * -1;
	double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z;


	if (obj.refl == DIFF) {
		// Lambertian

		// use QMC to sample the next direction
		//double r1 = 2.*PI*hal(d3 - 1, i), r2 = hal(d3 + 0, i);
		double r1 = 2.*PI*Random(), r2 = Random();
		double r2s = sqrt(r2);
		Vec w = nl, u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).norm();
		Vec v = w % u, d = (u*cos(r1)*r2s + v * sin(r1)*r2s + w * sqrt(1 - r2)).norm();

		if (m) {
			// eye ray
			// store the measurment point
			HPoint &hp = hitPoints[pixel];
			hp.used = true;
			hp.importance = f.mul(adj);
			hp.pos = x;
			hp.nrm = n;
			hp.pix = pixel;
			hitPoints[pixel] = hp;
			int lightindex = (int)(sizeof(sph) / sizeof(Sphere)) - 1;
			if ((dpt == 1 || specularbounce) && id == lightindex)
				directillum[hp.pix] = directillum[hp.pix] + adj.mul(sph[lightindex].e);
			else
				directillum[hp.pix] = directillum[hp.pix] + DirectIllumination(x, nl, hp.f, Vec(Random(), Random(), Random()));
		}
		else
		{
			// photon ray
			// find neighboring measurement points and accumulate flux via progressive density estimation
			Vec hh = (x - hpbbox.minPoint) * hashCellSize;
			int ix = abs(int(hh.x)), iy = abs(int(hh.y)), iz = abs(int(hh.z));
			// strictly speaking, we should use #pragma omp critical here.
			// it usually works without an artifact due to the fact that photons are 
			// rarely accumulated to the same measurement points at the same time (especially with QMC).
			// it is also significantly faster.
			if (dpt > 1) {
				std::vector<HPoint*> &hp = hashGrid[hash(ix, iy, iz)];
				for (HPoint* hitpoint : hp) {
					Vec v = hitpoint->pos - x;
					if ((hitpoint->nrm.dot(n) > 1e-3) && (v.dot(v) <= radius2[hitpoint->pix])) {
						// unlike N in the paper, hitpoint->n stores "N / ALPHA" to make it an integer value
						double g = (photonNums[hitpoint->pix] * ALPHA + ALPHA) / (photonNums[hitpoint->pix] * ALPHA + 1.0);
						//double g = (photonNums[hitpoint->pix] + ALPHA) / (photonNums[hitpoint->pix] + 1.0);
						radius2[hitpoint->pix] = radius2[hitpoint->pix] * g;
						photonNums[hitpoint->pix]++;
						//photonNums[hitpoint->pix] += ALPHA;
						flux[hitpoint->pix] = (flux[hitpoint->pix] + hitpoint->f.mul(fl)*(1. / PI))*g;
					}
				}
			}

			//if (hal(d3 + 1, i)<p) trace(Ray(x, d), dpt, m, f.mul(fl)*(1. / p), adj, i);
			if (Random() < p) trace(Ray(x, d), dpt, m, f.mul(fl)*(1. / p), adj, i);
		}

	}
	else if (obj.refl == SPEC) {
		// mirror
		trace(Ray(x, r.d - n * 2.0*n.dot(r.d)), dpt, m, f.mul(fl), f.mul(adj), i, true, pixel);

	}
	else {
		// glass
		Ray lr(x, r.d - n * 2.0*n.dot(r.d));
		bool into = (n.dot(nl) > 0.0);
		double nc = 1.0, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.d.dot(nl), cos2t;

		// total internal reflection
		if ((cos2t = 1 - nnt * nnt*(1 - ddn * ddn)) < 0) return trace(lr, dpt, m, fl, adj, i, true, pixel);

		Vec td = (r.d*nnt - n * ((into ? 1 : -1)*(ddn*nnt + sqrt(cos2t)))).norm();
		double a = nt - nc, b = nt + nc, R0 = a * a / (b*b), c = 1 - (into ? -ddn : td.dot(n));
		double Re = R0 + (1 - R0)*c*c*c*c*c, P = Re; Ray rr(x, td); Vec fa = f.mul(adj);
		if (m) {
			// eye ray (trace both rays)
			//trace(lr, dpt, m, fl, fa*Re, i, true, pixel);
			//trace(rr, dpt, m, fl, fa*(1.0 - Re), i, true, pixel);
			P = 0.5 * Re + 0.25;
			(Random() < P) ? trace(lr, dpt, m, fl, fa*Re / P, i, true, pixel) : trace(rr, dpt, m, fl, fa*(1.0 - Re) / (1 - P), i, true, pixel);
		}
		else {
			// photon ray (pick one via Russian roulette)
			//(hal(d3 - 1, i)<P) ? trace(lr, dpt, m, fl, fa, i) : trace(rr, dpt, m, fl, fa, i);
			(Random() < P) ? trace(lr, dpt, m, fl, fa, i) : trace(rr, dpt, m, fl, fa, i);
		}
	}
}
*/

void TraceEyePath(const Scene &scene, const Ray &ray, int maxDepth, long long pixel) {
	Ray r = ray;
	bool deltaBoundEvent = false;
	Vec importance(1.0, 1.0, 1.0);
	for (int i = 0; i < maxDepth; ++i) {
		double t;
		std::shared_ptr<Shape> hitObj;
		if (!scene.Intersect(r, t, hitObj)) return;
		std::cout << "Hit" << std::endl;
		Vec hit = r.o + r.d * t;
		Vec n = hitObj->GetNorm(hit);
		Vec nl = n.dot(r.d) < 0 ? n : n * -1;
		std::shared_ptr<BSDF> bsdf = hitObj->GetBSDF(n, nl);
		Vec wi;
		double pdf;
		if (bsdf->IsDelta()) {
			Vec f = bsdf->Sample_f(-1 * r.d, &wi, &pdf, Vec(Random(), Random(), Random()));
			importance = f * std::abs(wi.dot(n)) * importance / pdf;
			r.o = hit;
			r.d = wi;
		}
		else {
			//Vec f = bsdf->Sample_f(-1 * r.d, &wi, &pdf, Vec(Random(), Random(), Random()));
			//importance = f * std::abs(wi.dot(n)) * importance / pdf;

			HPoint &hp = hitPoints[pixel];
			hp.used = true;
			hp.importance = importance;
			hp.pos = hit;
			hp.nrm = n;
			hp.pix = pixel;
			hp.outDir = -1 * r.d;
			hitPoints[pixel] = hp;
			if ((i == 0 || deltaBoundEvent) && hitObj->IsLight())
				directillum[hp.pix] = directillum[hp.pix] + importance * hitObj->GetEmission();
			else
				directillum[hp.pix] = directillum[hp.pix] + 
					DirectIllumination(scene, hit, -1 * r.d, bsdf, nl, hp.importance, Vec(Random(), Random(), Random()));

			return;
		}
	}
}

void TracePhoton(const Scene &scene, const Ray &ray, Vec photonFlux, int maxDepth) {
	Ray r = ray;
	bool deltaBoundEvent = false;
	for (int i = 0; i < maxDepth; ++i) {
		double t;
		std::shared_ptr<Shape> hitObj;
		if (!scene.Intersect(r, t, hitObj)) return;
		Vec hit = r.o + r.d * t;
		Vec n = hitObj->GetNorm(hit);
		Vec nl = n.dot(r.d) < 0 ? n : n * -1;
		std::shared_ptr<BSDF> bsdf = hitObj->GetBSDF(n, nl);
		Vec wi;
		double pdf;
		Vec f = bsdf->Sample_f(-1 * r.d, &wi, &pdf, Vec(Random(), Random(), Random()));
		Vec estimation = f * std::abs(wi.dot(n)) / pdf;
		if (bsdf->IsDelta()) {
			photonFlux = photonFlux * estimation;
			r.o = hit;
			r.d = wi;
		}
		else {
			// photon ray
			// find neighboring measurement points and accumulate flux via progressive density estimation
			Vec hh = (hit - hpbbox.minPoint) * hashCellSize;
			int ix = std::abs(int(hh.x)), iy = std::abs(int(hh.y)), iz = std::abs(int(hh.z));
			// strictly speaking, we should use #pragma omp critical here.
			// it usually works without an artifact due to the fact that photons are 
			// rarely accumulated to the same measurement points at the same time (especially with QMC).
			// it is also significantly faster.
			if (i > 1) {
				std::vector<HPoint*> &hp = hashGrid[hash(ix, iy, iz)];
				for (HPoint* hitpoint : hp) {
					Vec v = hitpoint->pos - hit;
					if ((hitpoint->nrm.dot(n) > 1e-3) && (v.dot(v) <= radius2[hitpoint->pix])) {
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
		}
	}
}



int main(int argc, char *argv[]) {
	// samps * 1000 photon paths will be traced
	int w = 1024, h = 768, samps = (argc == 2) ? std::max(atoll(argv[1]) / render_stage_number, (long long)1) : render_stage_number;

	radius2.resize(w * h);
	photonNums.resize(w * h);
	flux.resize(w * h);
	hitPoints.resize(w * h);
	directillum.resize(w * h);

	std::shared_ptr<Scene> scene = std::shared_ptr<Scene>(new Scene);

	// trace eye rays and store measurement points
	Ray cam(Vec(50, 48, 295.6), Vec(0, -0.042612, -1).norm());
	Vec cx = Vec(w*.5135 / h), cy = (cx%cam.d).norm()*.5135, *c = new Vec[w*h], vw;

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
	std::shared_ptr<Sphere> light0 = std::shared_ptr<Sphere>(new Sphere(8.0, Vec(50, 81.6 - 16.5, 81.6), Vec(0.3, 0.3, 0.3) * 100, Vec(), DIFF));


	std::shared_ptr<Light> light = std::shared_ptr(new SphereLight(new))
	scene->AddShape(std::dynamic_pointer_cast<Shape>(light0));
	//Light
	scene->AddLight(std::dynamic_pointer_cast<Light>(light0));

	auto shapes = scene->GetShapes();
	auto lights = scene->GetLights();

	std::cout << (shapes[10] == std::dynamic_pointer_cast<Shape>(lights[0])) << std::endl;


	for (int render_stage = 0; render_stage < samps; ++render_stage) {
		for (int y = 0; y < h; y++) {
			fprintf(stderr, "\rHitPointPass %5.2f%%", 100.0*y / (h - 1));
#pragma omp parallel for schedule(dynamic, 1)
			for (int x = 0; x < w; x++) {
				int pixel = x + y * w;
				hitPoints[pixel].used = false;
				double u = Random() - 0.5, v = Random() - 0.5;
				Vec d = cx * ((x + 0.5 + u) / w - 0.5) + cy * (-(y + 0.5 + v) / h + 0.5) + cam.d;
				//trace(Ray(cam.o + d * 140, d.norm()), 0, true, Vec(), Vec(1, 1, 1), 0, false, pixel);
				Ray ray(cam.o + d * 140, d.norm());
				TraceEyePath(*scene, ray, 20, pixel);
			}
		}

		BuildHashGrid(w, h, render_stage);
		photonTotalNum = samps;
		vw = Vec(1, 1, 1);
		{
			fprintf(stderr, "\n");
			double p = 100.*(render_stage + 1) / photonTotalNum;
			fprintf(stderr, "\rPhotonPass %5.2f%%", p);
			int m = render_stage_number * render_stage;
			Ray ray;
			Vec photonFlux;
#pragma omp parallel for schedule(dynamic, 1)
			for (int j = 0; j < render_stage_number; j++) {
				GenratePhoton(*scene, &ray, &photonFlux);
				TracePhoton(*scene, ray, photonFlux, 20);
				//GenratePhoton(&r, &f, m + j);
				//trace(r, 0, 0 > 1, f, vw, m + j);
			}
			fprintf(stderr, "\n");
		}
		ClearHashGrid();
	}

	// density estimation
	/*
	for (int i = 0; i < flux.size(); ++i) {
		std::cout << directillum[i].x << " " << directillum[i].x << std::endl;
	}*/

	for (int i = 0; i < flux.size(); ++i) {
		c[i] = c[i] + flux[i] * (1.0 / (PI * radius2[i] * photonTotalNum * render_stage_number))
			+ directillum[i] / samps;
	}

	// save the image after tone mapping and gamma correction
	FILE* f = fopen("cornellbox.png", "wb");
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
