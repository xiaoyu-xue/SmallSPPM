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
#include "brute_force.h"
#include "kdtree_accel.h"
#include "sppm.h"
#include "diffuse.h"
#include "mirror.h"
#include "glass.h"
#include "constant_texture.h"

#include "transform.h"

#include "cornellbox.h"

#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#define _CRTDBG_MAP_ALLOC

//#define DEBUG_TRANSMIT

//const real ALPHA = 0.66666667;
const real ALPHA = 0.75;
const int64  render_stage_number = 200000;

void TestSppm(int argc, char* argv[]) {
	//clock_t begin = clock();

	//int w = 1024, h = 768;
	//int nIterations = (argc == 2) ? atol(argv[1]) : 256; //(argc == 2) ? std::max(atoll(argv[1]) / render_stage_number, (int64)1) : render_stage_number;

	//std::shared_ptr<Film> film = std::shared_ptr<Film>(new Film(w, h));

	//std::shared_ptr<Scene> scene = std::shared_ptr<Scene>(new Scene);

	//// trace eye rays and store measurement points
	////Ray cam(Vec3(50, 48, 295.6), Vec3(0, -0.042612, -1).norm());
	////Vec3 cx = Vec3(w*.5135 / h), cy = (cx%cam.d).norm()*.5135, *c = new Vec3[w*h], vw;

	//Vec3 camPos(50, 52, 295.6f);
	//Vec3 cz(0, -0.042612f, -1);
	////Vec3 cz(0, 0, -1);
	////real filmDis = cz.length();
	//real filmDis = cz.Length();
	//Vec3 cx = Vec3(w * .5135f / h, 0, 0).Norm();
	//Vec3 cy = (cx.Cross(cz)).Norm();
	//real fovy = 28.7993f;

	////std::cout << camPos + cz << std::endl << cy << std::endl;

	//std::shared_ptr<Camera> camera = std::shared_ptr<Camera>(new PinHoleCamera(film, camPos, cz, cx, cy, fovy, filmDis));
	//std::shared_ptr<Sampler> randomSampler = std::shared_ptr<Sampler>(new RandomSampler(123));
	//std::shared_ptr<Sampler> haltonSampler = std::shared_ptr<Sampler>(new HaltonSampler(w, h));
	//std::shared_ptr<Sampler> regularHaltonSampler = std::shared_ptr<Sampler>(new RegularHaltonSampler());;
	//std::shared_ptr<SamplerEnum> samplerEnum = std::shared_ptr<SamplerEnum>(new SamplerEnum());
	//std::shared_ptr<SamplerEnum> haltonSamplerEnum = std::shared_ptr<SamplerEnum>(new HaltonEnum((unsigned)w, (unsigned)h));
	//std::shared_ptr<Integrator> integrator =
	//	std::shared_ptr<Integrator>(new SPPM(nIterations, render_stage_number, 20, 1.0, ALPHA, false, haltonSampler, haltonSamplerEnum));
	//fprintf(stderr, "Load Scene ...\n");
	////scene->SetCamera(cam, cx, cy);
	//scene->SetCamera(camera);
	//scene->AddShape(std::shared_ptr<Shape>(new Sphere(1e5f, Vec3(1e5f + 1, 40.8f, 81.6f), Vec3(), Vec3(.75f, .25f, .25f), DIFF)));//Left
	//scene->AddShape(std::shared_ptr<Shape>(new Sphere(1e5f, Vec3(-1e5f + 99, 40.8f, 81.6f), Vec3(), Vec3(.25f, .25f, .75f), DIFF)));//Right
	//scene->AddShape(std::shared_ptr<Shape>(new Sphere(1e5f, Vec3(50, 40.8f, 1e5), Vec3(), Vec3(.75f, .75f, .75f), DIFF)));//Back
	////scene->AddShape(std::shared_ptr<Shape>(new Sphere(1e5, Vec3(50, 40.8, -1e5 + 170), Vec3(), Vec3(.75, .75, .75), DIFF)));//Frnt
	//scene->AddShape(std::shared_ptr<Shape>(new Sphere(1e5f, Vec3(50, 1e5, 81.6f), Vec3(), Vec3(.75f, .75f, .75f), DIFF)));//Botm
	//scene->AddShape(std::shared_ptr<Shape>(new Sphere(1e5f, Vec3(50, -1e5f + 81.6f, 81.6f), Vec3(), Vec3(.75f, .75f, .75f), DIFF)));//Top
	//scene->AddShape(std::shared_ptr<Shape>(new Sphere(16.5f, Vec3(27, 16.5f, 47), Vec3(), Vec3(1, 1, 1), REFR)));//Glass
	////**scene->AddShape(std::shared_ptr<Shape>(new Sphere(7.0, Vec3(27, 16.5, 47), Vec3(), Vec3(.25, .25, .75), DIFF)));//Mirr
	//scene->AddShape(std::shared_ptr<Shape>(new Sphere(16.5f, Vec3(73, 26.5f, 78), Vec3(), Vec3(1, 1, 1), REFR)));//Glass
	//scene->AddShape(std::shared_ptr<Shape>(new Sphere(9.5f, Vec3(53, 9.5f, 88), Vec3(), Vec3(1, 1, 1), REFR)));//Glass
	//scene->AddShape(std::shared_ptr<Shape>(new Sphere(9.5f, Vec3(23, 0.0f, 98), Vec3(), Vec3(1, 1, 1), DIFF)));//DIFF
	//std::shared_ptr<Shape> lightShape = std::shared_ptr<Shape>(new Sphere(8.f, Vec3(50, 81.6f - 16.5f, 81.6f), Vec3(0.3f, 0.3f, 0.3f) * 100, Vec3(), DIFF));//Lite
	//std::shared_ptr<Light> light0 = std::shared_ptr<Light>(new AreaLight(lightShape));
	//scene->AddLight(light0);
	//scene->Initialize();
	//film->SetFileName("cornellbox29.bmp");
	//std::shared_ptr<Renderer> renderer = std::shared_ptr<Renderer>(new Renderer(scene, integrator, film));
	//renderer->Render();
	//clock_t end = clock();
	//std::cout << "cost time: " << (end - begin) / 1000.0 / 60.0 << " min" << std::endl;
}

void TestSPPM2(int argc, char* argv[]) {
	clock_t begin = clock();

	int w = 1024, h = 768;
	int nIterations = (argc == 2) ? atol(argv[1]) : 256;

	std::shared_ptr<Film> film = std::shared_ptr<Film>(new Film(w, h));
	std::shared_ptr<Scene> scene = std::shared_ptr<Scene>(new Scene);
	Vec3 camPos(50, 52, 295.6f);
	Vec3 cz(0, -0.042612f, -1);

	real filmDis = cz.Length();
	Vec3 cx = Vec3(w * .5135f / h, 0, 0).Norm();
	Vec3 cy = (cx.Cross(cz)).Norm();
	real fovy = 28.7993f;

	std::shared_ptr<Camera> camera = std::shared_ptr<Camera>(new PinHoleCamera(film, camPos, cz, cx, cy, fovy, filmDis));
	std::shared_ptr<Sampler> randomSampler = std::shared_ptr<Sampler>(new RandomSampler(123));
	std::shared_ptr<Sampler> haltonSampler = std::shared_ptr<Sampler>(new HaltonSampler(w, h));
	std::shared_ptr<Sampler> regularHaltonSampler = std::shared_ptr<Sampler>(new RegularHaltonSampler());;
	std::shared_ptr<SamplerEnum> samplerEnum = std::shared_ptr<SamplerEnum>(new SamplerEnum());
	std::shared_ptr<SamplerEnum> haltonSamplerEnum = std::shared_ptr<SamplerEnum>(new HaltonEnum((unsigned)w, (unsigned)h));
	std::shared_ptr<Integrator> integrator =
		std::shared_ptr<Integrator>(new SPPM(nIterations, render_stage_number, 20, 1.0, ALPHA, false, haltonSampler, haltonSamplerEnum));
	fprintf(stderr, "Load Scene ...\n");

	std::shared_ptr<Accelerator> accelerator = std::shared_ptr<Accelerator>(new BruteForce());
	scene->SetCamera(camera);
	scene->SetAccelerator(accelerator);
	//texture
	std::shared_ptr<Texture<Vec3>> redConstant = std::shared_ptr<Texture<Vec3>>(new ConstantTexture<Vec3>(Vec3(.75f, .25f, .25f)));
	std::shared_ptr<Texture<Vec3>> blueConstant = std::shared_ptr<Texture<Vec3>>(new ConstantTexture<Vec3>(Vec3(.25f, .25f, .75f)));
	std::shared_ptr<Texture<Vec3>> whiteConstant = std::shared_ptr<Texture<Vec3>>(new ConstantTexture<Vec3>(Vec3(.75f, .75f, .75f)));
	std::shared_ptr<Texture<Vec3>> fullWhiteConstant = std::shared_ptr<Texture<Vec3>>(new ConstantTexture<Vec3>(Vec3(1, 1, 1)));


	//Left
	std::shared_ptr<Shape> leftWallShape = std::shared_ptr<Shape>(new Sphere(nullptr, nullptr, 1e5f, Vec3(1e5f + 1, 40.8f, 81.6f)));
	std::shared_ptr<Material> leftWallMaterial = std::shared_ptr<Material>(new DiffuseMaterial(redConstant));
	std::shared_ptr<Primitive> leftWall = std::shared_ptr<Primitive>(new GeometryPrimitive(leftWallShape, leftWallMaterial));
	scene->AddPrimitive(leftWall);

	//Right
	std::shared_ptr<Shape> rightWallShape = std::shared_ptr<Shape>(new Sphere(nullptr, nullptr, 1e5f, Vec3(-1e5f + 99, 40.8f, 81.6f)));
	std::shared_ptr<Material> rightWallMaterial = std::shared_ptr<Material>(new DiffuseMaterial(blueConstant));
	std::shared_ptr<Primitive> rightWall = std::shared_ptr<Primitive>(new GeometryPrimitive(rightWallShape, rightWallMaterial));
	scene->AddPrimitive(rightWall);

	//Back
	std::shared_ptr<Shape> backWallShape = std::shared_ptr<Shape>(new Sphere(nullptr, nullptr, 1e5f, Vec3(50, 40.8f, 1e5)));
	std::shared_ptr<Material> backWallMaterial = std::shared_ptr<Material>(new DiffuseMaterial(whiteConstant));
	std::shared_ptr<Primitive> backWall = std::shared_ptr<Primitive>(new GeometryPrimitive(backWallShape, backWallMaterial));
	scene->AddPrimitive(backWall);

	//Botom
	std::shared_ptr<Shape> botomShape = std::shared_ptr<Shape>(new Sphere(nullptr, nullptr, 1e5f, Vec3(50, 1e5, 81.6f)));
	std::shared_ptr<Material> botomMaterial = std::shared_ptr<Material>(new DiffuseMaterial(whiteConstant));
	std::shared_ptr<Primitive> botom = std::shared_ptr<Primitive>(new GeometryPrimitive(botomShape, botomMaterial));
	scene->AddPrimitive(botom);

	//Top
	std::shared_ptr<Shape> topShape = std::shared_ptr<Shape>(new Sphere(nullptr, nullptr, 1e5f, Vec3(50, -1e5f + 81.6f, 81.6f)));
	std::shared_ptr<Material> topMaterial = std::shared_ptr<Material>(new DiffuseMaterial(whiteConstant));
	std::shared_ptr<Primitive> top = std::shared_ptr<Primitive>(new GeometryPrimitive(topShape, topMaterial));
	scene->AddPrimitive(top);

	//Diffuse Ball1
	std::shared_ptr<Shape> diffuseBallShape1 = std::shared_ptr<Shape>(new Sphere(nullptr, nullptr, 9.5f, Vec3(23, 0.0f, 98)));
	std::shared_ptr<Material> diffuseBallMaterial = std::shared_ptr<Material>(new DiffuseMaterial(fullWhiteConstant));
	std::shared_ptr<Primitive> diffuseBall1 = std::shared_ptr<Primitive>(new GeometryPrimitive(diffuseBallShape1, diffuseBallMaterial));
	scene->AddPrimitive(diffuseBall1);

	//Glass Ball1
	std::shared_ptr<Shape> glassBallShape1 = std::shared_ptr<Shape>(new Sphere(nullptr, nullptr, 16.5f, Vec3(73, 26.5f, 78)));
	std::shared_ptr<Material> glassBallMaterial = std::shared_ptr<Material>(new GlassMaterial(fullWhiteConstant, fullWhiteConstant));
	std::shared_ptr<Primitive> glassBall1 = std::shared_ptr<Primitive>(new GeometryPrimitive(glassBallShape1, glassBallMaterial));
	scene->AddPrimitive(glassBall1);

	//Diffuse Ball2 new Sphere(7.0, Vec3(27, 16.5, 47), Vec3(), Vec3(.25, .25, .75)
	std::shared_ptr<Shape> diffuseBallShape2 = std::shared_ptr<Shape>(new Sphere(nullptr, nullptr, 7.0, Vec3(27, 16.5, 47)));
	std::shared_ptr<Material> diffuseBallMaterial2 = std::shared_ptr<Material>(new DiffuseMaterial(blueConstant));
	std::shared_ptr<Primitive> diffuseBall2 = std::shared_ptr<Primitive>(new GeometryPrimitive(diffuseBallShape2, diffuseBallMaterial2));
	scene->AddPrimitive(diffuseBall2);


	//Glass Ball2 new Sphere(16.5f, Vec3(27, 16.5f, 47)
	std::shared_ptr<Shape> glassBallShape2 = std::shared_ptr<Shape>(new Sphere(nullptr, nullptr, 16.5f, Vec3(27, 16.5f, 47)));
	std::shared_ptr<Primitive> glassBall2 = std::shared_ptr<Primitive>(new GeometryPrimitive(glassBallShape2, glassBallMaterial));
	scene->AddPrimitive(glassBall2);

	//Glass Ball3 new Sphere(9.5f, Vec3(53, 9.5f, 88)
	std::shared_ptr<Shape> glassBallShape3 = std::shared_ptr<Shape>(new Sphere(nullptr, nullptr, 9.5f, Vec3(53, 9.5f, 88)));
	std::shared_ptr<Primitive> glassBall3 = std::shared_ptr<Primitive>(new GeometryPrimitive(glassBallShape3, glassBallMaterial));
	scene->AddPrimitive(glassBall3);

	//Light
	std::shared_ptr<Texture<Vec3>> lightTexture = std::shared_ptr<Texture<Vec3>>(new ConstantTexture<Vec3>(Vec3(0, 0, 0)));
	std::shared_ptr<Shape> lightShape = std::shared_ptr<Shape>(new Sphere(nullptr, nullptr, 8.f, Vec3(50, 81.6f - 16.5f, 81.6f)));//Lite
	std::shared_ptr<Light> light0 = std::shared_ptr<Light>(new AreaLight(lightShape, Vec3(0.3f, 0.3f, 0.3f) * 100));
	std::shared_ptr<Material> lightMaterial = std::shared_ptr<Material>(new DiffuseMaterial(lightTexture));
	std::shared_ptr<Primitive> lightPrimitive = std::shared_ptr<Primitive>(new GeometryPrimitive(lightShape, lightMaterial, light0));
	scene->AddLight(lightPrimitive);

	scene->Initialize();
	film->SetFileName("cornellbox31.bmp");
	std::shared_ptr<Renderer> renderer = std::shared_ptr<Renderer>(new Renderer(scene, integrator, film));
	renderer->Render();
	clock_t end = clock();
	std::cout << "cost time: " << (end - begin) / 1000.0 / 60.0 << " min" << std::endl;
}

void TestSPPM3(int argc, char* argv[]) {
	clock_t begin = clock();

	int w = 1024, h = 768;
	int nIterations = (argc == 2) ? atol(argv[1]) : 256;

	std::shared_ptr<Film> film = std::shared_ptr<Film>(new Film(w, h));
	std::shared_ptr<Scene> scene = std::shared_ptr<Scene>(new Scene);
	Vec3 camPos(50, 52, 295.6f);
	Vec3 cz(0, -0.042612f, -1);

	real filmDis = cz.Length();
	Vec3 cx = Vec3(w * .5135f / h, 0, 0).Norm();
	Vec3 cy = (cx.Cross(cz)).Norm();
	real fovy = 28.7993f;

	std::shared_ptr<Camera> camera = std::shared_ptr<Camera>(new PinHoleCamera(film, camPos, cz, cx, cy, fovy, filmDis));
	std::shared_ptr<Sampler> randomSampler = std::shared_ptr<Sampler>(new RandomSampler(123));
	std::shared_ptr<Sampler> haltonSampler = std::shared_ptr<Sampler>(new HaltonSampler(w, h));
	std::shared_ptr<Sampler> regularHaltonSampler = std::shared_ptr<Sampler>(new RegularHaltonSampler());;
	std::shared_ptr<SamplerEnum> samplerEnum = std::shared_ptr<SamplerEnum>(new SamplerEnum());
	std::shared_ptr<SamplerEnum> haltonSamplerEnum = std::shared_ptr<SamplerEnum>(new HaltonEnum((unsigned)w, (unsigned)h));
	std::shared_ptr<Integrator> integrator =
		std::shared_ptr<Integrator>(new SPPM(nIterations, render_stage_number, 20, 1.0, ALPHA, false, haltonSampler, haltonSamplerEnum));
	fprintf(stderr, "Load Scene ...\n");

	CornellBoxSphere::SetScene(scene);
	//CornellBoxTriangle::SetScene(scene);

	std::shared_ptr<Accelerator> accelerator = std::shared_ptr<Accelerator>(new BruteForce());
	scene->SetCamera(camera);
	scene->SetAccelerator(accelerator);



	////texture
	//std::shared_ptr<Texture> redConstant = std::shared_ptr<Texture>(new ConstantTexture(Vec3(.75f, .25f, .25f)));
	//std::shared_ptr<Texture> blueConstant = std::shared_ptr<Texture>(new ConstantTexture(Vec3(.25f, .25f, .75f)));
	//std::shared_ptr<Texture> whiteConstant = std::shared_ptr<Texture>(new ConstantTexture(Vec3(.75f, .75f, .75f)));
	//std::shared_ptr<Texture> fullWhiteConstant = std::shared_ptr<Texture>(new ConstantTexture(Vec3(1, 1, 1)));


	////Left
	//Transform leftWallObjectToWorld = Transform::Translate(Vec3(1e5f + 1, 40.8f, 81.6f));
	//Transform leftWallWorldToObject = Inverse(leftWallObjectToWorld);
	//std::shared_ptr<Shape> leftWallShape = 
	//	std::shared_ptr<Shape>(new Sphere(&leftWallObjectToWorld, &leftWallWorldToObject, 1e5f, Vec3(0, 0, 0), true));
	//std::shared_ptr<Material> leftWallMaterial = std::shared_ptr<Material>(new DiffuseMaterial(redConstant));
	//std::shared_ptr<Primitive> leftWall = std::shared_ptr<Primitive>(new GeometryPrimitive(leftWallShape, leftWallMaterial));
	//scene->AddPrimitive(leftWall);

	////Right
	//Transform rightWallObjectToWorld = Transform::Translate(Vec3(-1e5f + 99, 40.8f, 81.6f));
	//Transform rightWallWorldToObject = Inverse(rightWallObjectToWorld);
	//std::shared_ptr<Shape> rightWallShape = 
	//	std::shared_ptr<Shape>(new Sphere(&rightWallObjectToWorld, &rightWallWorldToObject, 1e5f, Vec3(0, 0, 0), true));
	//std::shared_ptr<Material> rightWallMaterial = std::shared_ptr<Material>(new DiffuseMaterial(blueConstant));
	//std::shared_ptr<Primitive> rightWall = std::shared_ptr<Primitive>(new GeometryPrimitive(rightWallShape, rightWallMaterial));
	//scene->AddPrimitive(rightWall);

	////Back
	//Transform backWallObjectToWorld = Transform::Translate(Vec3(50, 40.8f, 1e5));
	//Transform backWallWorldToObject = Inverse(backWallObjectToWorld);
	//std::shared_ptr<Shape> backWallShape = 
	//	std::shared_ptr<Shape>(new Sphere(&backWallObjectToWorld, &backWallWorldToObject, 1e5f, Vec3(0, 0, 0), true));
	//std::shared_ptr<Material> backWallMaterial = std::shared_ptr<Material>(new DiffuseMaterial(whiteConstant));
	//std::shared_ptr<Primitive> backWall = std::shared_ptr<Primitive>(new GeometryPrimitive(backWallShape, backWallMaterial));
	//scene->AddPrimitive(backWall);

	////Botom
	//Transform botomWallObjectToWorld = Transform::Translate(Vec3(50, 1e5, 81.6f));
	//Transform botomWallWorldToObject = Inverse(botomWallObjectToWorld);
	//std::shared_ptr<Shape> botomShape = 
	//	std::shared_ptr<Shape>(new Sphere(&botomWallObjectToWorld, &botomWallWorldToObject, 1e5f, Vec3(0, 0, 0), true));
	//std::shared_ptr<Material> botomMaterial = std::shared_ptr<Material>(new DiffuseMaterial(whiteConstant));
	//std::shared_ptr<Primitive> botom = std::shared_ptr<Primitive>(new GeometryPrimitive(botomShape, botomMaterial));
	//scene->AddPrimitive(botom);

	////Top
	//Transform topWallObjectToWorld = Transform::Translate(Vec3(50, -1e5f + 81.6f, 81.6f));
	//Transform topWallWorldToObject = Inverse(topWallObjectToWorld);
	//std::shared_ptr<Shape> topShape = 
	//	std::shared_ptr<Shape>(new Sphere(&topWallObjectToWorld, &topWallWorldToObject, 1e5f, Vec3(0, 0, 0), true));
	//std::shared_ptr<Material> topMaterial = std::shared_ptr<Material>(new DiffuseMaterial(whiteConstant));
	//std::shared_ptr<Primitive> top = std::shared_ptr<Primitive>(new GeometryPrimitive(topShape, topMaterial));
	//scene->AddPrimitive(top);

	////Diffuse Ball1
	//Transform diffuseBall1ObjectToWorld = Transform::Translate(Vec3(23, 0.0f, 98));
	//Transform diffuseBall1WorldToObject = Inverse(diffuseBall1ObjectToWorld);
	//std::shared_ptr<Shape> diffuseBallShape1 = 
	//	std::shared_ptr<Shape>(new Sphere(&diffuseBall1ObjectToWorld, &diffuseBall1WorldToObject, 9.5f, Vec3(0, 0, 0)));
	//std::shared_ptr<Material> diffuseBallMaterial = std::shared_ptr<Material>(new DiffuseMaterial(fullWhiteConstant));
	//std::shared_ptr<Primitive> diffuseBall1 = std::shared_ptr<Primitive>(new GeometryPrimitive(diffuseBallShape1, diffuseBallMaterial));
	//scene->AddPrimitive(diffuseBall1);

	////Glass Ball1
	//Transform glassBall1ObjectToWorld = Transform::Translate(Vec3(73, 26.5f, 78));
	//Transform glassBall1WorldToObject = Inverse(glassBall1ObjectToWorld);
	//std::shared_ptr<Shape> glassBallShape1 = 
	//	std::shared_ptr<Shape>(new Sphere(&glassBall1ObjectToWorld, &glassBall1WorldToObject, 16.5f, Vec3(0, 0, 0)));
	//std::shared_ptr<Material> glassBallMaterial = std::shared_ptr<Material>(new GlassMaterial(fullWhiteConstant, fullWhiteConstant));
	//std::shared_ptr<Primitive> glassBall1 = std::shared_ptr<Primitive>(new GeometryPrimitive(glassBallShape1, glassBallMaterial));
	//scene->AddPrimitive(glassBall1);

	////Diffuse Ball2 new Sphere(7.0, Vec3(27, 16.5, 47), Vec3(), Vec3(.25, .25, .75)
	////Transform diffuseBall2ObjectToWorld = Transform::Translate(Vec3(27, 16.5, 47));
	////Transform diffuseBall2WorldToObject = Inverse(diffuseBall2ObjectToWorld);
	////std::shared_ptr<Shape> diffuseBallShape2 = 
	////	std::shared_ptr<Shape>(new Sphere(&diffuseBall2ObjectToWorld, &diffuseBall2WorldToObject, 7.0, Vec3(0, 0, 0)));
	////std::shared_ptr<Material> diffuseBallMaterial2 = std::shared_ptr<Material>(new DiffuseMaterial(blueConstant));
	////std::shared_ptr<Primitive> diffuseBall2 = std::shared_ptr<Primitive>(new GeometryPrimitive(diffuseBallShape2, diffuseBallMaterial2));
	////scene->AddPrimitive(diffuseBall2);


	////Glass Ball2 new Sphere(16.5f, Vec3(27, 16.5f, 47)
	//Transform glassBall2ObjectToWorld = Transform::Translate(Vec3(27, 16.5f, 47));
	//Transform glassBall2WorldToObject = Inverse(glassBall2ObjectToWorld);
	//std::shared_ptr<Shape> glassBallShape2 = 
	//	std::shared_ptr<Shape>(new Sphere(&glassBall2ObjectToWorld, &glassBall2WorldToObject, 16.5f, Vec3(0, 0, 0)));
	//std::shared_ptr<Primitive> glassBall2 = std::shared_ptr<Primitive>(new GeometryPrimitive(glassBallShape2, glassBallMaterial));
	//scene->AddPrimitive(glassBall2);

	////Glass Ball3 new Sphere(9.5f, Vec3(53, 9.5f, 88)
	//Transform glassBall3ObjectToWorld = Transform::Translate(Vec3(53, 9.5f, 88));
	//Transform glassBall3WorldToObject = Inverse(glassBall3ObjectToWorld);
	//std::shared_ptr<Shape> glassBallShape3 = 
	//	std::shared_ptr<Shape>(new Sphere(&glassBall3ObjectToWorld, &glassBall3WorldToObject, 9.5f, Vec3(0, 0, 0)));
	//std::shared_ptr<Primitive> glassBall3 = std::shared_ptr<Primitive>(new GeometryPrimitive(glassBallShape3, glassBallMaterial));
	//scene->AddPrimitive(glassBall3);

	////Light Sphere(8.f, Vec3(50, 81.6f - 16.5f, 81.6f)
	//Transform sphereLightObjectToWorld = Transform::Translate(Vec3(0, 0, 0));
	//Transform sphereLightWorldToObject = Inverse(sphereLightObjectToWorld);
	//std::shared_ptr<Texture> lightTexture = std::shared_ptr<Texture>(new ConstantTexture(Vec3(0, 0, 0)));
	//std::shared_ptr<Shape> lightShape = 
	//	std::shared_ptr<Shape>(new Sphere(&sphereLightObjectToWorld, &sphereLightWorldToObject, 8.f, Vec3(50, 81.6f - 16.5f, 81.6f)));//Lite
	//std::shared_ptr<Light> light0 = std::shared_ptr<Light>(new AreaLight(lightShape, Vec3(0.3f, 0.3f, 0.3f) * 100));
	//std::shared_ptr<Material> lightMaterial = std::shared_ptr<Material>(new DiffuseMaterial(lightTexture));
	//std::shared_ptr<Primitive> lightPrimitive = std::shared_ptr<Primitive>(new GeometryPrimitive(lightShape, lightMaterial, light0));
	//scene->AddLight(lightPrimitive);



	scene->Initialize();
	film->SetFileName("cornellbox50.bmp");
	std::shared_ptr<Renderer> renderer = std::shared_ptr<Renderer>(new Renderer(scene, integrator, film));
	renderer->Render();
	clock_t end = clock();
	std::cout << "cost time: " << (end - begin) / 1000.0 / 60.0 << " min" << std::endl;
}

void TestProjection() {
	int w = 1024, h = 768;
	std::shared_ptr<Film> film = std::shared_ptr<Film>(new Film(w, h));
	std::shared_ptr<Scene> scene = std::shared_ptr<Scene>(new Scene);

	Vec3 camPos(50, 52, 295.6f);
	//Vec3 cz(0, -0.042612f, -1);
	Vec3 cz(0, 0, -1);
	real filmDis = cz.Length();
	Vec3 cx = Vec3(w * .5135f / h, 0, 0).Norm();
	Vec3 cy = (cx.Cross(cz)).Norm();
	real fovy = 28.7993f;
	std::shared_ptr<ProjectiveCamera> camera = std::shared_ptr<ProjectiveCamera>(new PinHoleCamera(film, camPos, cz, cx, cy, fovy, filmDis));
	
	//raster to ndc
	{
		std::cout << "raster to ndc:" << std::endl;
		Vec3 rasterPoint(0, 0, 0.5);
		Vec3 rasterToNDC = camera->RasterToNDC(rasterPoint);
		std::cout << rasterToNDC << std::endl;
	}

	//ndc to camera
	{
		std::cout << "ndc to camera:" << std::endl;
		Vec3 cameraPoint = camera->WorldToCamera(camera->GetFilm()->LU);
		//cameraPoint.z = 1000;
		std::cout << cameraPoint << std::endl;
		Vec3 cameraToNDC = camera->CameraToNDC(cameraPoint);
		std::cout << cameraToNDC << std::endl;

		Vec3 ndcToCamera = camera->NDCToCamera(Vec3(-1, 1, -1));
		std::cout << ndcToCamera << std::endl;
	}

}

void TestSPPM4(int argc, char* argv[]) {
	clock_t begin = clock();

	int resX = 1024, resY = 1024;
	int nIterations = (argc == 2) ? atol(argv[1]) : 256;

	std::shared_ptr<Film> film = std::shared_ptr<Film>(new Film(resX, resY));
	std::shared_ptr<Scene> scene = std::shared_ptr<Scene>(new Scene);
	Vec3 camPos(0, 0, 3);
	Vec3 cz(0, 0, -1);
	Vec3 cx = Vec3(1, 0, 0);
	Vec3 cy = Vec3(0, 1, 0);
	real filmDis = 1;
	real fovy = 53.13010235415597f;

	std::shared_ptr<Camera> camera = std::shared_ptr<Camera>(new PinHoleCamera(film, camPos, cz, cx, cy, fovy, filmDis));
	std::shared_ptr<Sampler> randomSampler = std::shared_ptr<Sampler>(new RandomSampler(123));
	std::shared_ptr<Sampler> haltonSampler = std::shared_ptr<Sampler>(new HaltonSampler(resX, resY));
	std::shared_ptr<Sampler> regularHaltonSampler = std::shared_ptr<Sampler>(new RegularHaltonSampler());;
	std::shared_ptr<SamplerEnum> samplerEnum = std::shared_ptr<SamplerEnum>(new SamplerEnum());
	std::shared_ptr<SamplerEnum> haltonSamplerEnum = std::shared_ptr<SamplerEnum>(new HaltonEnum((unsigned)resX, (unsigned)resY));
	real alpha = 0.66666667;
	std::shared_ptr<Integrator> integrator =
		std::shared_ptr<Integrator>(new SPPM(nIterations, render_stage_number, 20, 0.05, alpha, false, haltonSampler, haltonSamplerEnum));
	fprintf(stderr, "Load Scene ...\n");

	CornellBoxTriangle2::SetScene(scene);

	//std::shared_ptr<Accelerator> accelerator = std::shared_ptr<Accelerator>(new BruteForce());
	std::shared_ptr<Accelerator> accelerator = std::shared_ptr<Accelerator>(new KdTreeAccel(scene->GetPrimitives()));
	scene->SetCamera(camera);
	scene->SetAccelerator(accelerator);

	scene->Initialize();
	film->SetFileName("cornellboxTexture6.bmp");
	std::shared_ptr<Renderer> renderer = std::shared_ptr<Renderer>(new Renderer(scene, integrator, film));
	renderer->Render();
	clock_t end = clock();
	std::cout << "cost time: " << (end - begin) / 1000.0 / 60.0 << " min" << std::endl;

}


void TestSPPM5(int argc, char* argv[]) {
	clock_t begin = clock();

	int resX = 1024, resY = 1024;
	int nIterations = (argc == 2) ? atol(argv[1]) : 256;

	std::shared_ptr<Film> film = std::shared_ptr<Film>(new Film(resX, resY));
	std::shared_ptr<Scene> scene = std::shared_ptr<Scene>(new Scene);
	Vec3 camPos(0, 0, 3);
	Vec3 cz(0, 0, -1);
	Vec3 cx = Vec3(1, 0, 0);
	Vec3 cy = Vec3(0, 1, 0);
	real filmDis = 1;
	real fovy = 53.13010235415597f;

	std::shared_ptr<Camera> camera = std::shared_ptr<Camera>(new PinHoleCamera(film, camPos, cz, cx, cy, fovy, filmDis));
	std::shared_ptr<Sampler> randomSampler = std::shared_ptr<Sampler>(new RandomSampler(123));
	std::shared_ptr<Sampler> haltonSampler = std::shared_ptr<Sampler>(new HaltonSampler(resX, resY));
	std::shared_ptr<Sampler> regularHaltonSampler = std::shared_ptr<Sampler>(new RegularHaltonSampler());;
	std::shared_ptr<SamplerEnum> samplerEnum = std::shared_ptr<SamplerEnum>(new SamplerEnum());
	std::shared_ptr<SamplerEnum> haltonSamplerEnum = std::shared_ptr<SamplerEnum>(new HaltonEnum((unsigned)resX, (unsigned)resY));
	real alpha = 0.66666667;
	std::shared_ptr<Integrator> integrator =
		std::shared_ptr<Integrator>(new SPPM(nIterations, render_stage_number, 50, 0.05, alpha, false, haltonSampler, haltonSamplerEnum));
	fprintf(stderr, "Load Scene ...\n");

	CornellBoxMesh::SetScene(scene);


	//std::shared_ptr<Accelerator> accelerator = std::shared_ptr<Accelerator>(new BruteForce());
	std::shared_ptr<Accelerator> accelerator = std::shared_ptr<Accelerator>(new KdTreeAccel(scene->GetPrimitives()));
	scene->SetCamera(camera);
	scene->SetAccelerator(accelerator);

	scene->Initialize();
	film->SetFileName("cornellboxMeshObj9.bmp");
	std::shared_ptr<Renderer> renderer = std::shared_ptr<Renderer>(new Renderer(scene, integrator, film));
	renderer->Render();
	clock_t end = clock();
	std::cout << "cost time: " << (end - begin) / 1000.0 / 60.0 << " min" << std::endl;

}

int main(int argc, char *argv[]) {
	//AABB aabb;
	//aabb = Union(Union(aabb, Vec3(-1, -2, -2)), Vec3(1, 2, 3));
	//std::cout << aabb << std::endl;

	//std::shared_ptr<Texture<Vec3>> imageTexture1 =
	//	std::shared_ptr<Texture<Vec3>>(new ImageTexture<Vec3>("..\\texture_images\\checkboard.bmp"));
	//imageTexture1->Sample(Vec2(0.1, 0.2));
	TestSPPM5(argc, argv);
	//TestSPPM3(argc, argv);
	//_CrtDumpMemoryLeaks();
}
