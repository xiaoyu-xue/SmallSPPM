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
const real ALPHA = 0.75;
const int64  render_stage_number = 20000;




int main(int argc, char *argv[]) {


	clock_t begin = clock();

	int w = 1024, h = 768;
	int nIterations = (argc == 2) ? atol(argv[1]) : 256; //(argc == 2) ? std::max(atoll(argv[1]) / render_stage_number, (int64)1) : render_stage_number;

	std::shared_ptr<Film> film = std::shared_ptr<Film>(new Film(w, h));

	std::shared_ptr<Scene> scene = std::shared_ptr<Scene>(new Scene);

	// trace eye rays and store measurement points
	//Ray cam(Vec3(50, 48, 295.6), Vec3(0, -0.042612, -1).norm());
	//Vec3 cx = Vec3(w*.5135 / h), cy = (cx%cam.d).norm()*.5135, *c = new Vec3[w*h], vw;

	Vec3 camPos(50, 52, 295.6f), cz(0, -0.042612f, -1);
	//real filmDis = cz.length();
	real filmDis = 1.0f;
	Vec3 cx = Vec3(w * .5135f / h, 0, 0).Norm();
	Vec3 cy = (cx.Cross(cz)).Norm();
	real fovy = 28.7993f;

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
	scene->AddShape(std::shared_ptr<Shape>(new Sphere(1e5f, Vec3(1e5f + 1, 40.8f, 81.6f), Vec3(), Vec3(.75f, .25f, .25f), DIFF)));//Left
	scene->AddShape(std::shared_ptr<Shape>(new Sphere(1e5f, Vec3(-1e5f + 99, 40.8f, 81.6f), Vec3(), Vec3(.25f, .25f, .75f), DIFF)));//Right
	scene->AddShape(std::shared_ptr<Shape>(new Sphere(1e5f, Vec3(50, 40.8f, 1e5), Vec3(), Vec3(.75f, .75f, .75f), DIFF)));//Back
	//scene->AddShape(std::shared_ptr<Shape>(new Sphere(1e5, Vec3(50, 40.8, -1e5 + 170), Vec3(), Vec3(.75, .75, .75), DIFF)));//Frnt
	scene->AddShape(std::shared_ptr<Shape>(new Sphere(1e5f, Vec3(50, 1e5, 81.6f), Vec3(), Vec3(.75f, .75f, .75f), DIFF)));//Botm
	scene->AddShape(std::shared_ptr<Shape>(new Sphere(1e5f, Vec3(50, -1e5f + 81.6f, 81.6f), Vec3(), Vec3(.75f, .75f, .75f), DIFF)));//Top
	//scene->AddShape(std::shared_ptr<Shape>(new Sphere(16.5f, Vec3(27, 16.5f, 47), Vec3(), Vec3(1, 1, 1), DIFF)));//Glass
	//**scene->AddShape(std::shared_ptr<Shape>(new Sphere(7.0, Vec3(27, 16.5, 47), Vec3(), Vec3(.25, .25, .75), DIFF)));//Mirr
	//scene->AddShape(std::shared_ptr<Shape>(new Sphere(16.5f, Vec3(73, 26.5f, 78), Vec3(), Vec3(1, 1, 1), DIFF)));//Glass
	//scene->AddShape(std::shared_ptr<Shape>(new Sphere(9.5f, Vec3(53, 9.5f, 88), Vec3(), Vec3(1, 1, 1), DIFF)));//Glass
	//scene->AddShape(std::shared_ptr<Shape>(new Sphere(9.5f, Vec3(23, 0.0f, 98), Vec3(), Vec3(1, 1, 1), DIFF)));//DIFF
	std::shared_ptr<Shape> lightShape = std::shared_ptr<Shape>(new Sphere(8.f, Vec3(50, 81.6f - 16.5f, 81.6f), Vec3(0.3f, 0.3f, 0.3f) * 100, Vec3(), DIFF));//Lite
	std::shared_ptr<Light> light0 = std::shared_ptr<Light>(new AreaLight(lightShape));
	scene->AddLight(light0);
	scene->Initialize();
	film->SetFileName("cornellbox11.bmp");
	std::shared_ptr<Renderer> renderer = std::shared_ptr<Renderer>(new Renderer(scene, integrator, film));
	renderer->Render();
	
	clock_t end = clock();


	std::cout << "cost time: " << (end - begin) / 1000.0 / 60.0 << " min" << std::endl;


	//_CrtDumpMemoryLeaks();
}
