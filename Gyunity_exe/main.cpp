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
#include "common/Platform.h"
#include "common/Core.h"
#include "sampler\Sampler.h"
#include "sampler\Halton.h"
#include "sampler\Sobol.h"
#include "sampler\SamplerEnum.h"
#include "sampler\HaltonEnum.h"
#include "sampler\SobolEnum.h"
#include "math\Linagl.h"
#include "core\Scene.h"
#include "core\Film.h"
#include "shape\sphere.h"
#include "integrator\Integrator.h"
#include "core\Renderer.h"
#include "light\Light.h"
#include "light\AreaLight.h"
#include "camera\PinholeCamera.h"
#include "accelerator\BruteForceAccel.h"
#include "accelerator\KdTreeAccel.h"
#include "accelerator\BVHAccel.h"
#include "integrator\SPPM.h"
#include "material\Diffuse.h"
#include "material\Mirror.h"
#include "material\Glass.h"
#include "texture\ConstantTexture.h"
#include "math\Transform.h"
#include "integrator\PT.h"
#include "math\Filter.h"
#include "integrator\VPT.h"
#include "CornellBox.h"
#include "integrator\LT.h"
#include "common/Config.h"
#include <windows.h>
#include <cstdlib>
#include "integrator/BDPT.h"

using namespace Gyunity;

#define _CRTDBG_MAP_ALLOC

//#define DEBUG_TRANSMIT

//const real ALPHA = 0.66666667;
const real ALPHA = 0.75;
const int64  render_stage_number = 300000;

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
	//std::cout << "\ncost time: " << (end - begin) / 1000.0 / 60.0 << " min" << std::endl;
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
	scene->AddLight(light0);
	scene->AddPrimitive(lightPrimitive);

	scene->Initialize();
	film->SetFileName("cornellbox31.bmp");
	std::shared_ptr<Renderer> renderer = std::shared_ptr<Renderer>(new Renderer(scene, camera, integrator, film));
	renderer->Render();
	clock_t end = clock();
	std::cout << "\ncost time: " << (end - begin) / 1000.0 / 60.0 << " min" << std::endl;
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

	//std::shared_ptr<Accelerator> accelerator = std::shared_ptr<Accelerator>(new BruteForce());
	std::shared_ptr<Accelerator> accelerator = std::shared_ptr<Accelerator>(new BVHAccel(scene->GetPrimitives()));

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
	std::shared_ptr<Renderer> renderer = std::shared_ptr<Renderer>(new Renderer(scene, camera, integrator, film));
	renderer->Render();
	clock_t end = clock();
	std::cout << "\ncost time: " << (end - begin) / 1000.0 / 60.0 << " min" << std::endl;
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
		Vec3 cameraPoint = camera->WorldToCamera(camera->GetFilm()->mLU);
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

	int resX = 1024 * 2, resY = 1024 * 2;
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
		std::shared_ptr<Integrator>(new SPPM(nIterations, render_stage_number, 20, 0.05, alpha, true, haltonSampler, haltonSamplerEnum));
	fprintf(stderr, "Load Scene ...\n");

	CornellBoxTriangle2::SetScene(scene);

	//std::shared_ptr<Accelerator> accelerator = std::shared_ptr<Accelerator>(new BruteForce());
	std::shared_ptr<Accelerator> accelerator = std::shared_ptr<Accelerator>(new KdTreeAccel(scene->GetPrimitives()));
	scene->SetAccelerator(accelerator);

	scene->Initialize();
	film->SetFileName("cornellboxTexture6.bmp");
	std::shared_ptr<Renderer> renderer = std::shared_ptr<Renderer>(new Renderer(scene, camera, integrator, film));
	renderer->Render();
	clock_t end = clock();
	std::cout << "\ncost time: " << (end - begin) / 1000.0 / 60.0 << " min" << std::endl;

}


void TestSPPM5(int argc, char* argv[]) {


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
		std::shared_ptr<Integrator>(new SPPM(nIterations, render_stage_number, 20, 0.05, alpha, false, haltonSampler, haltonSamplerEnum, true));
	fprintf(stderr, "Load Scene ...\n");

	CornellBoxMesh::SetScene(scene);
	//CornellBoxWater::SetScene(scene);
	//CornellBoxTriangle2::SetScene(scene);
	//EnvironmentMapScene::SetScene(scene);
	//CornellBoxHeartSurface::SetScene(scene);
	//HeartSurfaceEnvironmentMapScene::SetScene(scene);

	//std::shared_ptr<Accelerator> accelerator = std::shared_ptr<Accelerator>(new BruteForce());
	//std::shared_ptr<Accelerator> accelerator = std::shared_ptr<Accelerator>(new KdTreeAccel(scene->GetPrimitives()));
	std::shared_ptr<Accelerator> accelerator = std::shared_ptr<Accelerator>(new BVHAccel(scene->GetPrimitives()));
	scene->SetAccelerator(accelerator);

	scene->Initialize();
	film->SetFileName("CornellBoxSppmTest.bmp");
	std::shared_ptr<Renderer> renderer = std::shared_ptr<Renderer>(new Renderer(scene, camera, integrator, film));
	clock_t begin = clock();
	renderer->Render();
	clock_t end = clock();
	std::cout << "\ncost time: " << (end - begin) / 1000.0 / 60.0 << " min" << std::endl;

}

void TestHashGrid() {
	Rng rng;
	std::vector<Vec3> points;
	for (int i = 0; i < 12345; ++i) {
		points.push_back(Vec3(rng.GetFloat(), rng.GetFloat(), rng.GetFloat()));
	}
	real searchRadius = 0.15;
	HashGrid<int> hashGrid;
	for (int i = 0; i < points.size(); ++i) {
		hashGrid.AddPoint(std::move(std::pair<Vec3, int>(points[i], i)), 0.15f);
	}
	hashGrid.BuildHashGrid(searchRadius + Eps);
	Vec3 testPoint = Vec3(0.12f, 0.34f, 0.56f);


	{
		int naiveNum = 0;
		for (int i = 0; i < points.size(); ++i) {
			//std::cout << points[i] << std::endl;
			if (Distance(points[i], testPoint) < searchRadius) {
				naiveNum++;
				//std::cout << points[i] << std::endl;
			}
		}
		std::cout << naiveNum << std::endl;
	}
	{
		int hashNum = 0;
		std::vector<int> grid = hashGrid.GetGrid(testPoint);
		for (int i = 0; i < grid.size(); ++i) {
			if (Distance(points[grid[i]], testPoint) < searchRadius) {
				hashNum++;
				//std::cout << points[grid[i]] << std::endl;
			}
		}
		std::cout << hashNum << std::endl;
	}


}

void TestSPPM6(int argc, char* argv[]) {


	int resX = 1024, resY = 1024;
	int nIterations = (argc == 2) ? atol(argv[1]) : 256;

	std::shared_ptr<Film> film = std::shared_ptr<Film>(new Film(resX, resY, new BoxFilter()));
	std::shared_ptr<Scene> scene = std::shared_ptr<Scene>(new Scene);
	Vec3 camPos(0, 0, 3);
	Vec3 lookAt(0, 0, 0);
	Vec3 up(0, 1, 0);
	real filmDis = 1;
	real fovy = 53.13010235415597f;

	std::shared_ptr<Camera> camera = std::shared_ptr<Camera>(new PinHoleCamera(film, camPos, lookAt, up, fovy, filmDis));
	std::shared_ptr<Sampler> randomSampler = std::shared_ptr<Sampler>(new RandomSampler(123));
	std::shared_ptr<Sampler> haltonSampler = std::shared_ptr<Sampler>(new HaltonSampler(resX, resY));
	std::shared_ptr<Sampler> regularHaltonSampler = std::shared_ptr<Sampler>(new RegularHaltonSampler());;
	std::shared_ptr<SamplerEnum> samplerEnum = std::shared_ptr<SamplerEnum>(new SamplerEnum());
	std::shared_ptr<SamplerEnum> haltonSamplerEnum = std::shared_ptr<SamplerEnum>(new HaltonEnum((unsigned)resX, (unsigned)resY));
	real alpha = 0.66666667;
	std::shared_ptr<Integrator> integrator =
		std::shared_ptr<Integrator>(new SPPM(nIterations, render_stage_number, 20, 0.05, alpha, false, haltonSampler, haltonSamplerEnum, true));

	GYT_Print("Load Scene ...\n");

	CornellBoxMesh::SetScene(scene);
	//CornellBoxWater::SetScene(scene);
	//CornellBoxTriangle2::SetScene(scene);
	//EnvironmentMapScene::SetScene(scene);
	//CornellBoxHeartSurface::SetScene(scene);
	//HeartSurfaceEnvironmentMapScene::SetScene(scene);

	std::shared_ptr<Accelerator> accelerator = std::shared_ptr<Accelerator>(new BVHAccel(scene->GetPrimitives()));
	scene->SetAccelerator(accelerator);

	scene->Initialize();
	film->SetFileName("CornellBoxSppmTest3.bmp");
	std::shared_ptr<Renderer> renderer = std::shared_ptr<Renderer>(new Renderer(scene, camera, integrator, film));
	clock_t begin = clock();
	renderer->Render();
	clock_t end = clock();
	std::cout << "\ncost time: " << (end - begin) / 1000.0 / 60.0 << " min" << std::endl;

}

void TestVolPathTracing(int argc, char* argv[]) {

	int resX = 1024, resY = 1024;

	std::shared_ptr<Film> film = std::shared_ptr<Film>(new Film(resX, resY, new BoxFilter()));
	std::shared_ptr<Scene> scene = std::shared_ptr<Scene>(new Scene);
	Vec3 camPos(0, 0, 3);
	Vec3 cz(0, 0, -1);
	Vec3 cx = Vec3(1, 0, 0);
	Vec3 cy = Vec3(0, 1, 0);
	real filmDis = 1;
	real fovy = 53.13010235415597f;

	std::shared_ptr<Camera> camera = std::shared_ptr<Camera>(new PinHoleCamera(film, camPos, cz, cx, cy, fovy, filmDis));
	std::shared_ptr<Sampler> randomSampler = std::shared_ptr<Sampler>(new RandomSampler(123));
	std::shared_ptr<SamplerEnum> samplerEnum = std::shared_ptr<SamplerEnum>(new SamplerEnum());

	std::shared_ptr<Sampler> sobolSampler = std::shared_ptr<Sampler>(new SobolSampler());
	std::shared_ptr<SamplerEnum> sobolSamplerEnum = std::shared_ptr<SamplerEnum>(new SobolEnum(resX, resY));

	//std::shared_ptr<Integrator> integrator = std::shared_ptr<Integrator>(new VolPathTracing(10, 10, sobolSampler, sobolSamplerEnum));
	std::shared_ptr<Integrator> integrator = std::shared_ptr<Integrator>(new VolPathTracing(2000, 8, randomSampler, samplerEnum, false));

	fprintf(stderr, "Load Scene ...\n");

	//CornellBoxMesh::SetScene(scene);
	//CornellBoxTriangle2::SetScene(scene);
	//EnvironmentMapScene::SetScene(scene);
	//CornellBoxHeartSurface::SetScene(scene);
	//HeartSurfaceEnvironmentMapScene::SetScene(scene);
	CornellBoxMedium::SetScene(scene);

	std::shared_ptr<Accelerator> accelerator = std::shared_ptr<Accelerator>(new BVHAccel(scene->GetPrimitives()));
	scene->SetAccelerator(accelerator);

	scene->Initialize();
	film->SetFileName("MediaTestFogEquiAngular12.bmp");
	std::shared_ptr<Renderer> renderer = std::shared_ptr<Renderer>(new Renderer(scene, camera, integrator, film));
	clock_t begin = clock();
	renderer->Render();
	clock_t end = clock();
	std::cout << "\ncost time: " << (end - begin) / 1000.0 / 60.0 << " min" << std::endl;
}

void TestTransmittance() {
	std::shared_ptr<Scene> scene = std::shared_ptr<Scene>(new Scene);
	CornellBoxMedium::SetScene(scene);
	Intersection p0(Vec3(0, -0.99999, 0), Vec3(), MediumInterface());
	p0.mNormal = p0.mAbsNormal = p0.mGeometryNormal = Vec3(0, 1, 0);
	Intersection p1(Vec3(0, 0.965, 0), Vec3(), MediumInterface());
	p1.mNormal = p1.mAbsNormal = p1.mGeometryNormal = Vec3(0, -1, 0);
	VisibilityTester visibilityTester(p0, p1);
	std::shared_ptr<Accelerator> accelerator = std::shared_ptr<Accelerator>(new BVHAccel(scene->GetPrimitives()));
	scene->SetAccelerator(accelerator);
	scene->Initialize();

	std::shared_ptr<Sampler> randomSampler = std::shared_ptr<Sampler>(new RandomSampler(123));
	std::shared_ptr<SamplerEnum> samplerEnum = std::shared_ptr<SamplerEnum>(new SamplerEnum());
	RandomStateSequence rand(randomSampler, 123);
	Vec3 Tr = visibilityTester.Tr(*scene, rand);
	std::cout << Tr << std::endl;
	std::cout << std::exp(-0.6) << std::endl;
}


void TestPathTracing2(int argc, char* argv[]) {

	int resX = 1024, resY = 1024;

	std::shared_ptr<Film> film = std::shared_ptr<Film>(new Film(resX, resY, new BoxFilter()));
	std::shared_ptr<Scene> scene = std::shared_ptr<Scene>(new Scene);
	Vec3 camPos(0, 0, 3);
	Vec3 lookAt(0, 0, 0);
	Vec3 up(0, 1, 0);
	real filmDis = 1;
	real fovy = 53.13010235415597f;

	std::shared_ptr<Camera> camera = std::shared_ptr<Camera>(new PinHoleCamera(film, camPos, lookAt, up, fovy, filmDis));
	std::shared_ptr<Sampler> randomSampler = std::shared_ptr<Sampler>(new RandomSampler(123));
	std::shared_ptr<SamplerEnum> samplerEnum = std::shared_ptr<SamplerEnum>(new SamplerEnum());

	std::shared_ptr<Sampler> sobolSampler = std::shared_ptr<Sampler>(new SobolSampler());
	std::shared_ptr<SamplerEnum> sobolSamplerEnum = std::shared_ptr<SamplerEnum>(new SobolEnum(resX, resY));

	//std::shared_ptr<Integrator> integrator = std::shared_ptr<Integrator>(new PathTracing(100, 20, sobolSampler, sobolSamplerEnum));
	std::shared_ptr<Integrator> integrator = std::shared_ptr<Integrator>(new PathTracing(100, 1, randomSampler, samplerEnum));

	//fprintf(stderr, "Load Scene ...\n");
	GYT_Print("Load Scene ...\n");
	//CornellBoxMesh::SetScene(scene);
	//CornellBoxMeshEmpty::SetScene(scene);
	CornellBoxMeshWithSlope::SetScene(scene);
	//CornellBoxDiffuseBunny::SetScene(scene);
	//CornellBoxTriangle2::SetScene(scene);
	//EnvironmentMapScene::SetScene(scene);
	//CornellBoxHeartSurface::SetScene(scene);
	//HeartSurfaceEnvironmentMapScene::SetScene(scene);

	std::shared_ptr<Accelerator> accelerator = std::shared_ptr<Accelerator>(new BVHAccel(scene->GetPrimitives()));
	scene->SetAccelerator(accelerator);

	scene->Initialize();
	film->SetFileName("CornellBoxPT10.bmp");
	std::shared_ptr<Renderer> renderer = std::shared_ptr<Renderer>(new Renderer(scene, camera, integrator, film));
	clock_t begin = clock();
	renderer->Render();
	clock_t end = clock();

	GYT_Print("\nCost time: {:.6f} min\n", (end - begin) / 1000.0 / 60.0);
}

void RenderAdScene(int idx, real theta) {

	int resX = 2160, resY = 2160;
	std::shared_ptr<Film> film = std::shared_ptr<Film>(new Film(resX, resY, new BoxFilter()));
	std::shared_ptr<Scene> scene = std::shared_ptr<Scene>(new Scene);
	Vec3 camPos(0, 0, 3);
	Vec3 lookAt(0, 0, 0);
	Vec3 up(0, 1, 0);
	real filmDis = 1;
	real fovy = 53.13010235415597f;

	std::shared_ptr<Camera> camera = std::shared_ptr<Camera>(new PinHoleCamera(film, camPos, lookAt, up, fovy, filmDis));
	std::shared_ptr<Sampler> randomSampler = std::shared_ptr<Sampler>(new RandomSampler(123));
	std::shared_ptr<SamplerEnum> samplerEnum = std::shared_ptr<SamplerEnum>(new SamplerEnum());

	std::shared_ptr<Integrator> integrator = std::shared_ptr<Integrator>(new PathTracing(256, 10, randomSampler, samplerEnum));

	GYT_Print("Load Scene ...\n");

	//AdBoardScene::SetScene(scene, theta);
	AdBoardSceneV2 buildScene;
	buildScene.SetScene(scene, theta);


	std::shared_ptr<Accelerator> accelerator = std::shared_ptr<Accelerator>(new BVHAccel(scene->GetPrimitives()));
	scene->SetAccelerator(accelerator);

	scene->Initialize();
	film->SetFileName("./Result/" + std::to_string(idx) + ".png");
	std::shared_ptr<Renderer> renderer = std::shared_ptr<Renderer>(new Renderer(scene, camera, integrator, film));
	renderer->Render();

}

void TestLightTracing(int argc, char* argv[]) {


	int resX = 1024, resY = 1024;
	int nIterations = (argc == 2) ? atol(argv[1]) : 256;

	std::shared_ptr<Film> film = std::shared_ptr<Film>(new Film(resX, resY, new BoxFilter()));
	std::shared_ptr<Scene> scene = std::shared_ptr<Scene>(new Scene);
	Vec3 camPos(0, 0, 3);
	Vec3 lookAt(0, 0, 0);
	Vec3 up(0, 1, 0);
	real filmDis = 1;
	real fovy = 53.13010235415597f;

	std::shared_ptr<Camera> camera = std::shared_ptr<Camera>(new PinHoleCamera(film, camPos, lookAt, up, fovy, filmDis));
	std::shared_ptr<Sampler> randomSampler = std::shared_ptr<Sampler>(new RandomSampler(123));


	std::shared_ptr<Integrator> integrator = std::shared_ptr<Integrator>(new LightTracing(randomSampler, 15, 512));

	std::shared_ptr<Accelerator> accelerator = std::shared_ptr<Accelerator>(new BVHAccel(scene->GetPrimitives()));
	scene->SetAccelerator(accelerator);

	scene->Initialize();
	film->SetFileName("./Result/LT2.bmp");
	std::shared_ptr<Renderer> renderer = std::shared_ptr<Renderer>(new Renderer(scene, camera, integrator, film));
	clock_t begin = clock();
	renderer->Render();

	clock_t end = clock();
	std::cout << "\ncost time: " << (end - begin) / 1000.0 / 60.0 << " min" << std::endl;

}

void TestPathTracing(int argc, char* argv[]) {

	int resX = 1024, resY = 1024;
	int nIterations = (argc == 2) ? atol(argv[1]) : 256;

	std::shared_ptr<Film> film = std::shared_ptr<Film>(new Film(resX, resY, new BoxFilter()));
	std::shared_ptr<Scene> scene = std::shared_ptr<Scene>(new Scene);
	Vec3 camPos(0, 0, 3);
	Vec3 lookAt(0, 0, 0);
	Vec3 up(0, 1, 0);
	real filmDis = 1;
	real fovy = 53.13010235415597f;
	std::shared_ptr<Camera> camera = std::shared_ptr<Camera>(new PinHoleCamera(film, camPos, lookAt, up, fovy, filmDis));

	std::shared_ptr<Sampler> randomSampler = std::shared_ptr<Sampler>(new RandomSampler(123));
	std::shared_ptr<SamplerEnum> samplerEnum = std::shared_ptr<SamplerEnum>(new SamplerEnum());

	std::shared_ptr<Integrator> integrator = std::make_shared<PathTracing>(512, 15, randomSampler, samplerEnum);

	GYT_Print("Load Scene ...\n");
	SimpleCornellBox::SetScene(scene);

	std::shared_ptr<Accelerator> accelerator = std::shared_ptr<Accelerator>(new BVHAccel(scene->GetPrimitives()));
	scene->SetAccelerator(accelerator);

	scene->Initialize();
	film->SetFileName("./Result/PT.bmp");
	std::shared_ptr<Renderer> renderer = std::shared_ptr<Renderer>(new Renderer(scene, camera, integrator, film));
	clock_t begin = clock();
	renderer->Render();

	clock_t end = clock();
	std::cout << "\ncost time: " << (end - begin) / 1000.0 / 60.0 << " min" << std::endl;
}



void TestBDPT(int argc, char* argv[]) {

	int resX = 1024, resY = 1024;
	int nIterations = (argc == 2) ? atol(argv[1]) : 256;

	std::shared_ptr<Film> film = std::shared_ptr<Film>(new Film(resX, resY, new BoxFilter()));
	std::shared_ptr<Scene> scene = std::shared_ptr<Scene>(new Scene);
	Vec3 camPos(0, 0, 3);
	Vec3 lookAt(0, 0, 0);
	Vec3 up(0, 1, 0);
	real filmDis = 1;
	real fovy = 53.13010235415597f;
	std::shared_ptr<Camera> camera = std::shared_ptr<Camera>(new PinHoleCamera(film, camPos, lookAt, up, fovy, filmDis));

	std::shared_ptr<Sampler> randomSampler = std::shared_ptr<Sampler>(new RandomSampler(123));
	std::shared_ptr<SamplerEnum> samplerEnum = std::shared_ptr<SamplerEnum>(new SamplerEnum());

	std::shared_ptr<Integrator> integrator = std::make_shared<BDPT>(randomSampler, samplerEnum, 16, 128, false, false);

	GYT_Print("Load Scene ...\n");
	SimpleCornellBox::SetScene(scene);

	std::shared_ptr<Accelerator> accelerator = std::shared_ptr<Accelerator>(new BVHAccel(scene->GetPrimitives()));
	scene->SetAccelerator(accelerator);

	scene->Initialize();
	film->SetFileName("./Result/BDPT5.bmp");
	std::shared_ptr<Renderer> renderer = std::shared_ptr<Renderer>(new Renderer(scene, camera, integrator, film));
	clock_t begin = clock();
	renderer->Render();

	clock_t end = clock();
	std::cout << "\ncost time: " << (end - begin) / 1000.0 / 60.0 << " min" << std::endl;

}


int main(int argc, char *argv[]) {

	std::cout << GGXDistribution::RoughnessToAlpha(0.118) << std::endl;


	TestBDPT(argc, argv);
	//TestPathTracing(argc, argv);
	//TestLightTracing(argc, argv);
	//TestPathTracing2(argc, argv);
	//TestSPPM6(argc, argv);
	//TestPathTracing2(argc, argv);
	//TestPhathTracingargc, argv);
	//TestPathTracing2(argc, argv);
	//TestVolPathTracing(argc, argv);
	//TestTransmittance();
	

	//clock_t begin = clock();
	//for (int i = 0; i <= 360; i++) {
	//	GYT_Print("\nFrame: {}\n", i);
	//	RenderAdScene(i, (real)i);
	//}
	//clock_t end = clock();
	//std::cout << "\ncost time: " << (end - begin) / 1000.0 / 60.0 << " min" << std::endl;

 
	//RenderAdScene(0, 0);


	//RenderAdScene(281, (real)281);
	_CrtDumpMemoryLeaks();

	return 0;
}
