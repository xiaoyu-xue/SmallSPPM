// expanded smallppm (code is exactly the same as smallppm.cpp but with more comments)
#include "svpng.inc"
#include <math.h>   // smallppm, Progressive Photon Mapping by T. Hachisuka
#include <stdlib.h> // originally smallpt, a path tracer by Kevin Beason, 2008
#include <stdio.h>  // Usage: ./smallppm 100000 && xv image.ppm
#include <random>
#include <vector>
#include <crtdbg.h>
#include <iostream>
#include <algorithm>
#define _CRTDBG_MAP_ALLOC

#define PI ((double)3.14159265358979) // ^^^^^^:number of photons emitted
#define INV_PI ((double)0.31830988618379067154)
#define ALPHA ((double)0.7) // the alpha parameter of PPM

const int render_stage_number = 2000000;
const double PiOver2 = 1.57079632679489661923;
const double PiOver4 = 0.78539816339744830961;
const double eps = 10e-6;

std::mt19937_64 mRng(1234);
std::uniform_real_distribution<double> mDistdouble;

double random() {
	return mDistdouble(mRng);
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
	inline Vec operator/(double b) const { if (b == 0) return Vec(); else return Vec(x / b, y / b, z / b); }
	inline Vec mul(const Vec &b) const { return Vec(x * b.x, y * b.y, z * b.z); }
	inline Vec norm() { return (*this) * (1.0 / sqrt(x*x + y*y + z*z)); }
	inline double dot(const Vec &b) const { return x * b.x + y * b.y + z * b.z; }
	Vec operator%(Vec&b) const { return Vec(y*b.z - z*b.y, z*b.x - x*b.z, x*b.y - y*b.x); }
};

Vec operator*(double a, Vec b) { return Vec(a * b.x, a * b.y, a * b.z); }



void CoordinateSystem(const Vec &v1, Vec *v2, Vec *v3) {
	if (std::abs(v1.x) > std::abs(v1.y))
		*v2 = Vec(-v1.z, 0, v1.x) / std::sqrt(v1.x * v1.x + v1.z * v1.z);
	else
		*v2 = Vec(0, v1.z, -v1.y) / std::sqrt(v1.y * v1.y + v1.z * v1.z);
	*v3 = v1 % (*v2);
}

#define MAX(x, y) ((x > y) ? x : y)
struct AABB {
	Vec min, max; // axis aligned bounding box
	inline void fit(const Vec &p)
	{
		if (p.x<min.x)min.x = p.x; // min
		if (p.y<min.y)min.y = p.y; // min
		if (p.z<min.z)min.z = p.z; // min
		max.x = MAX(p.x, max.x);
		max.y = MAX(p.y, max.y);
		max.z = MAX(p.z, max.z);
	}
	inline void reset() {
		min = Vec(1e20, 1e20, 1e20);
		max = Vec(-1e20, -1e20, -1e20);
	}
};

struct HPoint {
	Vec f, pos, nrm, flux;
	double r2;
	unsigned int n; // n = N / ALPHA in the paper
	int pix;
	bool used;
};

std::vector<double> radius2;
std::vector<unsigned int> photon_num;
//std::vector<double> photon_num;
std::vector<Vec> flux;
std::vector<Vec> directillum;
std::vector<HPoint> HitPoints;


unsigned int num_hash, pixel_index, num_photon;
double hash_s; 
AABB hpbbox;
std::vector<std::vector<HPoint*>> hash_grid;

void clear_hash_grid() {
	for (auto &e : hash_grid) {
		std::vector<HPoint*>().swap(e);
	}
}

// spatial hash function
inline unsigned int hash(const int ix, const int iy, const int iz) {
	return (unsigned int)((ix * 73856093) ^ (iy * 19349663) ^ (iz * 83492791)) % num_hash;
}


void build_hash_grid(const int w, const int h, int stage_num) {
	// find the bounding box of all the measurement points
	hpbbox.reset();
	for (HPoint &hp : HitPoints) {
		if(hp.used) hpbbox.fit(hp.pos);
	}

	// heuristic for initial radius
	Vec ssize = hpbbox.max - hpbbox.min;
	double irad = ((ssize.x + ssize.y + ssize.z) / 3.0) / ((w + h) / 2.0) * 2.0;
	// determine hash table size
	// we now find the bounding box of all the measurement points inflated by the initial radius
	hpbbox.reset();
	int vphoton = 0;

	for (HPoint &hp : HitPoints) {
		if (!hp.used) continue;
		if (stage_num == 0) {
			radius2[hp.pix] = irad * irad;
			photon_num[hp.pix] = 0;
			flux[hp.pix] = Vec();
		}
		vphoton++;
		hpbbox.fit(hp.pos - irad);
		hpbbox.fit(hp.pos + irad);
	}

	// make each grid cell two times larger than the initial radius
	hash_s = 1.0 / (irad*2.0);
	num_hash = vphoton;

	// build the hash table

	hash_grid.resize(num_hash);

	for(HPoint &hp : HitPoints) {
		if (!hp.used) continue;
		Vec BMin = ((hp.pos - irad) - hpbbox.min) * hash_s;
		Vec BMax = ((hp.pos + irad) - hpbbox.min) * hash_s;
		for (int iz = abs(int(BMin.z)); iz <= abs(int(BMax.z)); iz++)
		{
			for (int iy = abs(int(BMin.y)); iy <= abs(int(BMax.y)); iy++)
			{
				for (int ix = abs(int(BMin.x)); ix <= abs(int(BMax.x)); ix++)
				{
					int hv = hash(ix, iy, iz);
					hash_grid[hv].push_back(&hp);
				}
			}
		}
	}
}


struct Ray { Vec o, d; Ray() {}; Ray(Vec o_, Vec d_) : o(o_), d(d_) {} };

enum Refl_t { DIFF, SPEC, REFR };  // material types, used in radiance()
struct Sphere {
	double rad; Vec p, e, c; Refl_t refl;
	Sphere(double r_, Vec p_, Vec e_, Vec c_, Refl_t re_) : rad(r_), p(p_), e(e_), c(c_), refl(re_) {}
	inline double intersect(const Ray &r) const {
		// ray-sphere intersection returns distance
		Vec op = p - r.o;
		double t, b = op.dot(r.d), det = b*b - op.dot(op) + rad*rad;
		if (det < 0) {
			return 1e20;
		}
		else {
			det = sqrt(det);
		}
		return (t = b - det) > 1e-4 ? t : ((t = b + det)>1e-4 ? t : 1e20);
	}
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


// tone mapping and gamma correction
int toInt(double x) {
	return int(pow(1 - exp(-x), 1 / 2.2) * 255 + .5);
}

// find the closet interection
inline bool intersect(const Ray &r, double &t, int &id) {
	int n = sizeof(sph) / sizeof(Sphere);
	double d, inf = 1e20; t = inf;
	for (int i = 0; i<n; i++) {
		d = sph[i].intersect(r);
		if (d<t) {
			t = d;
			id = i;
		}
	}
	return t<inf;
}

Vec ConcentricSampleDisk(const Vec &u) {
	// Map uniform random numbers to $[-1,1]^2$
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

// generate a photon ray from the point light source with QMC
void genp(Ray* pr, Vec* f, int i) {
	*f = Vec(2500, 2500, 2500) * 1.5 *(PI*4.0); // flux
	double p = 2.*PI*hal(0, i), t = 2.*acos(sqrt(1. - hal(1, i)));
	double st = sin(t);
	pr->d = Vec(cos(p)*st, cos(t), sin(p)*st);
	pr->o = Vec(50, 60, 85);
}

void genphoton(Ray *pr, Vec *f, int i) {
	int lightindex = (int)(sizeof(sph) / sizeof(Sphere))  - 1;
	const Sphere &light = sph[lightindex];
	Vec Le = light.e;
	//sample a position
	Vec pos = UniformSampleSphere(Vec(random(),random(),random())) * light.rad + light.p;
	double Pdfpos = 1.f / (4 * PI * light.rad * light.rad);
	Vec lightNorm = (pos - light.p).norm();
	Vec ss, ts;
	CoordinateSystem(lightNorm, &ss, &ts);
	Vec dirLocal = CosineSampleHemisphere(Vec(random(), random(), random()));
	double CosTheta = dirLocal.z;
	Vec lightDir = (ss * dirLocal.x + ts * dirLocal.y + lightNorm * dirLocal.z).norm();
	double Pdfdir = CosineHemispherePdf(CosTheta);
	pr->o = pos + lightDir * eps;
	pr->d = lightDir;
	*f = Le * CosTheta / (Pdfpos * Pdfdir);
}


Vec DirectIllumination(Vec hitpoint, Vec nl, Vec f, Vec u) {
	Vec e;
	int numSpheres = (int)(sizeof(sph) / sizeof(Sphere));
	for (int i = 0; i< numSpheres; i++) {
		const Sphere &s = sph[i];
		if (s.e.x <= 0 && s.e.y <= 0 && s.e.z <= 0) continue; // skip non-lights

		Vec sw = s.p - hitpoint, su = ((fabs(sw.x)>.1 ? Vec(0, 1) : Vec(1)) % sw).norm(), sv = sw%su;
		double cos_a_max = sqrt(1 - s.rad*s.rad / (hitpoint - s.p).dot(hitpoint - s.p));
		double eps1 = u.x, eps2 = u.y;
		double cos_a = 1 - eps1 + eps1*cos_a_max;
		double sin_a = sqrt(1 - cos_a*cos_a);
		double phi = 2 * PI*eps2;
		Vec l = su*cos(phi)*sin_a + sv*sin(phi)*sin_a + sw*cos_a;
		l = l.norm();
		double t;
		int id;
		if (intersect(Ray(hitpoint, l), t, id) && id == i) {  // shadow ray
			double omega = 2 * PI*(1 - cos_a_max);
			e = e + f.mul(s.e*l.dot(nl)*omega) * INV_PI;  // 1/pi for brdf
		}
	}
	return e;
}

void trace(const Ray &r, int dpt, bool m, const Vec &fl, const Vec &adj, int i, bool specularbounce = false, int pixel = -1)
{
	double t;
	int id;

	dpt++;
	if (!intersect(r, t, id) || (dpt >= 20))return;

	int d3 = dpt * 3;
	const Sphere &obj = sph[id];
	Vec x = r.o + r.d*t, n = (x - obj.p).norm(), f = obj.c;
	Vec nl = n.dot(r.d)<0 ? n : n*-1;
	double p = f.x>f.y&&f.x>f.z ? f.x : f.y>f.z ? f.y : f.z;


	if (obj.refl == DIFF) {
		// Lambertian

		// use QMC to sample the next direction
		//double r1 = 2.*PI*hal(d3 - 1, i), r2 = hal(d3 + 0, i);
		double r1 = 2.*PI*random(), r2 = random();
		double r2s = sqrt(r2);
		Vec w = nl, u = ((fabs(w.x)>.1 ? Vec(0, 1) : Vec(1)) % w).norm();
		Vec v = w%u, d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1 - r2)).norm();

		if (m) {
			// eye ray
			// store the measurment point
			HPoint &hp = HitPoints[pixel];
			hp.used = true;
			hp.f = f.mul(adj);
			hp.pos = x;
			hp.nrm = n;
			hp.pix = pixel;
			HitPoints[pixel] = hp;
			int lightindex = (int)(sizeof(sph) / sizeof(Sphere)) - 1;
			if ((dpt == 1 || specularbounce) && id == lightindex) 
				directillum[hp.pix] = directillum[hp.pix] + adj.mul(sph[lightindex].e);
			else 
				directillum[hp.pix] = directillum[hp.pix] + DirectIllumination(x, nl, hp.f, Vec(random(), random(), random()));
		}
		else
		{
			// photon ray
			// find neighboring measurement points and accumulate flux via progressive density estimation
			Vec hh = (x - hpbbox.min) * hash_s;
			int ix = abs(int(hh.x)), iy = abs(int(hh.y)), iz = abs(int(hh.z));
			// strictly speaking, we should use #pragma omp critical here.
			// it usually works without an artifact due to the fact that photons are 
			// rarely accumulated to the same measurement points at the same time (especially with QMC).
			// it is also significantly faster.
			if(dpt > 1){
				std::vector<HPoint*> &hp = hash_grid[hash(ix, iy, iz)];
				for(HPoint* hitpoint : hp){
					Vec v = hitpoint->pos - x;
					if ((hitpoint->nrm.dot(n) > 1e-3) && (v.dot(v) <= radius2[hitpoint->pix])) {
						// unlike N in the paper, hitpoint->n stores "N / ALPHA" to make it an integer value
						double g = (photon_num[hitpoint->pix]*ALPHA + ALPHA) / (photon_num[hitpoint->pix]*ALPHA + 1.0);
						//double g = (photon_num[hitpoint->pix] + ALPHA) / (photon_num[hitpoint->pix] + 1.0);
						radius2[hitpoint->pix] = radius2[hitpoint->pix]*g;
						photon_num[hitpoint->pix]++;
						//photon_num[hitpoint->pix] += ALPHA;
						flux[hitpoint->pix] = (flux[hitpoint->pix] + hitpoint->f.mul(fl)*(1. / PI))*g;
					}
				}
			}

			//if (hal(d3 + 1, i)<p) trace(Ray(x, d), dpt, m, f.mul(fl)*(1. / p), adj, i);
			if (random()<p) trace(Ray(x, d), dpt, m, f.mul(fl)*(1. / p), adj, i);
		}

	}
	else if (obj.refl == SPEC) {
		// mirror
		trace(Ray(x, r.d - n*2.0*n.dot(r.d)), dpt, m, f.mul(fl), f.mul(adj), i, true, pixel);

	}
	else {
		// glass
		Ray lr(x, r.d - n*2.0*n.dot(r.d));
		bool into = (n.dot(nl)>0.0);
		double nc = 1.0, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.d.dot(nl), cos2t;

		// total internal reflection
		if ((cos2t = 1 - nnt*nnt*(1 - ddn*ddn))<0) return trace(lr, dpt, m, fl, adj, i, true, pixel);

		Vec td = (r.d*nnt - n*((into ? 1 : -1)*(ddn*nnt + sqrt(cos2t)))).norm();
		double a = nt - nc, b = nt + nc, R0 = a*a / (b*b), c = 1 - (into ? -ddn : td.dot(n));
		double Re = R0 + (1 - R0)*c*c*c*c*c, P = Re; Ray rr(x, td); Vec fa = f.mul(adj);
		if (m) {
			// eye ray (trace both rays)
			//trace(lr, dpt, m, fl, fa*Re, i, true, pixel);
			//trace(rr, dpt, m, fl, fa*(1.0 - Re), i, true, pixel);
			P = 0.5 * Re + 0.25;
			(random() < P) ? trace(lr, dpt, m, fl, fa*Re/P, i, true, pixel) : trace(rr, dpt, m, fl, fa*(1.0 - Re)/(1-P), i, true, pixel);
		}
		else {
			// photon ray (pick one via Russian roulette)
			//(hal(d3 - 1, i)<P) ? trace(lr, dpt, m, fl, fa, i) : trace(rr, dpt, m, fl, fa, i);
			(random()<P) ? trace(lr, dpt, m, fl, fa, i) : trace(rr, dpt, m, fl, fa, i);
		}
	}
}





int main(int argc, char *argv[]) {
	// samps * 1000 photon paths will be traced
	int w = 1024, h = 768, samps = (argc == 2) ? MAX(atoll(argv[1]) / render_stage_number, 1) : render_stage_number;

	radius2.resize(w * h);
	photon_num.resize(w * h);
	flux.resize(w * h);
	HitPoints.resize(w * h);
	directillum.resize(w * h);

	// trace eye rays and store measurement points
	Ray cam(Vec(50, 48, 295.6), Vec(0, -0.042612, -1).norm());
	Vec cx = Vec(w*.5135 / h), cy = (cx%cam.d).norm()*.5135, *c = new Vec[w*h], vw;

	
	for (int render_stage = 0; render_stage < samps; ++render_stage) {
		for (int y = 0; y<h; y++) {
			fprintf(stderr, "\rHitPointPass %5.2f%%", 100.0*y / (h - 1));
#pragma omp parallel for schedule(dynamic, 1)
			for (int x = 0; x<w; x++) {
				int pixel = x + y * w;
				HitPoints[pixel].used = false;
				double u = random() - 0.5, v = random() - 0.5;
				Vec d = cx * ((x + 0.5 + u) / w - 0.5) + cy * (-(y + 0.5 + v) / h + 0.5) + cam.d;
				trace(Ray(cam.o + d * 140, d.norm()), 0, true, Vec(), Vec(1, 1, 1), 0, false, pixel);
			}
		}
		
		build_hash_grid(w, h, render_stage);
		num_photon = samps;
		vw = Vec(1, 1, 1);
		{
			fprintf(stderr, "\n");
			double p = 100.*(render_stage + 1) / num_photon;
			fprintf(stderr, "\rPhotonPass %5.2f%%", p);
			int m = render_stage_number * render_stage;
			Ray r;
			Vec f;
#pragma omp parallel for schedule(dynamic, 1)
			for (int j = 0; j < render_stage_number; j++) {
				//genp(&r, &f, m + j);
				genphoton(&r, &f, m + j);
				trace(r, 0, 0 > 1, f, vw, m + j);
			}
			fprintf(stderr, "\n");
		}
		clear_hash_grid();
	}
	
	// density estimation

	for (int i = 0; i < flux.size(); ++i) {
		c[i] = c[i] + flux[i] * (1.0 / (PI*radius2[i]*num_photon*render_stage_number)) 
			+ directillum[i] / samps;
	}

	// save the image after tone mapping and gamma correction
	FILE* f = fopen("image13.png", "wb"); 
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
