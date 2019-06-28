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

#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#define _CRTDBG_MAP_ALLOC

//#define DEBUG_TRANSMIT

//const real ALPHA = 0.66666667;
const real ALPHA = 0.70;
const int64  render_stage_number = 2000000;




int main(int argc, char *argv[]) {


	clock_t begin = clock();

	int w = 1024, h = 768;
	int nIterations = (argc == 2) ? atol(argv[1]) : 256; //(argc == 2) ? std::max(atoll(argv[1]) / render_stage_number, (int64)1) : render_stage_number;

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
	std::shared_ptr<Sampler> haltonSampler = std::shared_ptr<Sampler>(new HaltonSampler(w, h));
	std::shared_ptr<Sampler> regularHaltonSampler = std::shared_ptr<Sampler>(new RegularHaltonSampler());;
	std::shared_ptr<SamplerEnum> samplerEnum = std::shared_ptr<SamplerEnum>(new SamplerEnum());
	std::shared_ptr<SamplerEnum> haltonSamplerEnum = std::shared_ptr<SamplerEnum>(new HaltonEnum((unsigned)w, (unsigned)h));
	std::shared_ptr<Integrator> integrator = 
		std::shared_ptr<Integrator>(new SPPM(nIterations, render_stage_number, 20, 1.0, ALPHA, false, haltonSampler, haltonSamplerEnum));
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
	film->SetFileName("cornellbox7.bmp");
	std::shared_ptr<Renderer> renderer = std::shared_ptr<Renderer>(new Renderer(scene, integrator, film));
	renderer->Render();
	
	clock_t end = clock();


	std::cout << "cost time: " << (end - begin) / 1000.0 / 60.0 << " min" << std::endl;


	//_CrtDumpMemoryLeaks();
}
