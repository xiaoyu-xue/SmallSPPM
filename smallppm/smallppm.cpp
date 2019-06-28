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
#include "def.h"
#include "sampler.h"
#include "halton.h"
#include "sampler_enum.h"
#include "linagl.h"
#include "scene.h"
#include "sphere.h"
#include "integrator.h"
#include "renderer.h"
#include "light.h"
#include "arealight.h"
#include "phinhole.h"
#include "sppm.h"

#define _CRTDBG_MAP_ALLOC

//#define DEBUG_TRANSMIT

const real ALPHA = 0.66666667;
const int64  render_stage_number = 70000;

//class Shape;
//class BSDF;
//class Scene;
/*
std::mt19937_64 rng(1234);
std::uniform_real_distribution<real> uniform;
real Random() {
	return uniform(rng);
}
*/

/*
uint32 rand_int() noexcept {
	static uint32 x = 123456789, y = 362436069, z = 521288629, w = 88675123;
	uint32 t = x ^ (x << 11);
	x = y;
	y = z;
	z = w;
	return (w = (w ^ (w >> 19)) ^ (t ^ (t >> 8)));
}

real Rand() noexcept {
	 return rand_int() * (1.0f / 4294967296.0f);
 }*/





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

real hal(const int b, int j) {
	const int p = primes[b];
	real h = 0.0, f = 1.0 / (real)p, fct = f;
	while (j > 0) {
		h += rev(j % p, p) * fct; j /= p; fct *= f;
	}
	return h;
}*/





































//std::vector<real> radius2;
//std::vector<int64> photonNums;
//std::vector<Vec> flux;
//std::vector<Vec> directillum;
//std::vector<HPoint> hitPoints;
//std::vector<std::vector<HPoint*>> hashGrid;
//std::vector<Spinlock> hashGridSpinlocks;
//int64 hashNum, pixelIndex;
//real hashCellSize;
//AABB hpbbox;





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

/*
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
*/









int main(int argc, char *argv[]) {
	
	clock_t begin = clock();

	int w = 1024, h = 768;
	int nIterations = (argc == 2) ? atol(argv[1]) : 256; //(argc == 2) ? std::max(atoll(argv[1]) / render_stage_number, (int64)1) : render_stage_number;

	//std::cout << nIterations << std::endl;

	//Initialize(w, h, 0.8);
	//hashGridSpinlocks.resize(w * h);

	std::shared_ptr<Film> film = std::shared_ptr<Film>(new Film(w, h));

	std::shared_ptr<Scene> scene = std::shared_ptr<Scene>(new Scene);

	// trace eye rays and store measurement points
	//Ray cam(Vec(50, 48, 295.6), Vec(0, -0.042612, -1).norm());
	//Vec cx = Vec(w*.5135 / h), cy = (cx%cam.d).norm()*.5135, *c = new Vec[w*h], vw;

	Vec camPos(50, 52, 295.6), cz(0, -0.042612, -1);
	//real filmDis = cz.length();
	real filmDis = 1.0;
	Vec cx = Vec(w * .5135 / h).norm();
	Vec cy = (cx % cz).norm();
	real fovy = 28.7993;

	//std::cout << camPos + cz << std::endl << cy << std::endl;
	
	std::shared_ptr<Camera> camera = std::shared_ptr<Camera>(new PinHoleCamera(film, camPos, cz, cx, cy, fovy, filmDis));
	std::shared_ptr<Sampler> randomSampler = std::shared_ptr<Sampler>(new RandomSampler(123));
	std::shared_ptr<Sampler> haltonSampler = std::shared_ptr<Sampler>(new HaltonSampler());
	std::shared_ptr<SamplerEnum> haltonSamplerEnum = std::shared_ptr<SamplerEnum>(new HaltonEnum((unsigned)w, (unsigned)h));
	std::shared_ptr<Integrator> integrator = std::shared_ptr<Integrator>(new SPPM(nIterations, render_stage_number, 20, 0.8, ALPHA, true, haltonSampler, haltonSamplerEnum));
	fprintf(stderr, "Load Scene ...\n");
	//scene->SetCamera(cam, cx, cy);
	scene->SetCamera(camera);
	scene->AddShape(std::shared_ptr<Shape>(new Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(), Vec(.75, .25, .25), DIFF)));//Left
	scene->AddShape(std::shared_ptr<Shape>(new Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(), Vec(.25, .25, .75), DIFF)));//Right
	scene->AddShape(std::shared_ptr<Shape>(new Sphere(1e5, Vec(50, 40.8, 1e5), Vec(), Vec(.75, .75, .75), DIFF)));//Back
	//scene->AddShape(std::shared_ptr<Shape>(new Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Vec(), Vec(.75, .75, .75), DIFF)));//Frnt
	scene->AddShape(std::shared_ptr<Shape>(new Sphere(1e5, Vec(50, 1e5, 81.6), Vec(), Vec(.75, .75, .75), DIFF)));//Botm
	scene->AddShape(std::shared_ptr<Shape>(new Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(), Vec(.75, .75, .75), DIFF)));//Top
	scene->AddShape(std::shared_ptr<Shape>(new Sphere(16.5, Vec(27, 16.5, 47), Vec(), Vec(1, 1, 1), REFR)));//Mirr
	//scene->AddShape(std::shared_ptr<Shape>(new Sphere(7.0, Vec(27, 16.5, 47), Vec(), Vec(.25, .25, .75), DIFF)));//Mirr
	scene->AddShape(std::shared_ptr<Shape>(new Sphere(16.5, Vec(73, 26.5, 78), Vec(), Vec(1, 1, 1), REFR)));//Glass
	scene->AddShape(std::shared_ptr<Shape>(new Sphere(9.5, Vec(53, 9.5, 88), Vec(), Vec(1, 1, 1), REFR)));//Glass
	scene->AddShape(std::shared_ptr<Shape>(new Sphere(9.5, Vec(23, 0.0, 98), Vec(), Vec(1, 1, 1), DIFF)));//Glass
	std::shared_ptr<Shape> lightShape = std::shared_ptr<Shape>(new Sphere(8.0, Vec(50, 81.6 - 16.5, 81.6), Vec(0.3, 0.3, 0.3) * 100, Vec(), DIFF));//Lite
	std::shared_ptr<Light> light0 = std::shared_ptr<Light>(new AreaLight(lightShape));
	scene->AddLight(light0);
	scene->Initialize();
	film->SetFileName("cornellbox.bmp");
	std::shared_ptr<Renderer> renderer = std::shared_ptr<Renderer>(new Renderer(scene, integrator, film));
	renderer->Render();
	//film->ConvertBmpToPng("cornellbox31.bmp");

	Vector3 vec1(1, 2, 3), vec2(9, 8, 7);
	std::cout << vec1 << std::endl;
	std::cout << vec2 << std::endl;
	std::cout << vec1.Dot(vec2) << " " << Dot(vec1, vec2) << std::endl;
	std::cout << vec1.Cross(vec2) << " " << Cross(vec1, vec2) << std::endl;
	std::cout << vec1 * vec2 << std::endl;
	std::cout << vec1.Normal() << std::endl;
	std::cout << vec1 / 3.0 << std::endl;
	std::cout << Distance(vec1, vec2) << std::endl;
	Vector3 vec3;
	vec3 = vec1;
	std::cout << vec3 << std::endl;
	clock_t end = clock();

	std::cout << "cost time: "<< (end - begin) / 1000.0 / 60.0 <<" min"<< std::endl;



	//test
	/*
	Intersection isect;
	isect.n = Vec(0, 1, 0);
	isect.nl = isect.n * -1;
	isect.wo = Vec(1, -1, 0) * -1;
	Vec wi;
	real pdf;
	TransmissionBSDF bsdf(isect);
	bsdf.Sample_f(isect.wo, &wi, &pdf, Vec(0.5, 0.5, 0.5));
	std::cout << "------------------------" << std::endl;

	bool into = isect.n.dot(isect.nl) > 0;                // Ray from outside going in?
	real nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = (isect.wo * -1).dot(isect.nl), cos2t;
	cos2t = 1 - nnt * nnt*(1 - ddn * ddn);
	Vec tdir = ((isect.wo * -1)*nnt - isect.n * ((into ? 1 : -1)*(ddn*nnt + sqrt(cos2t)))).norm();
	real a = nt - nc, b = nt + nc, R0 = a * a / (b*b), c = 1 - (into ? -ddn : tdir.dot(isect.n));
	real Re = R0 + (1 - R0)*c*c*c*c*c;
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
