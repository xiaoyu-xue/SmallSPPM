#pragma once
#include "def.h"
#include "scene.h"
#include "sphere.h"
#include "diffuse.h"
#include "mirror.h"
#include "glass.h"
#include "arealight.h"
#include "constant_texture.h"
#include "triangle.h"

NAMESPACE_BEGIN

class CornellBoxSphere {
public:
	static void SetScene(std::shared_ptr<Scene>& scene) {
		//texture
		std::shared_ptr<Texture> redConstant = std::shared_ptr<Texture>(new ConstantTexture(Vec3(.75f, .25f, .25f)));
		std::shared_ptr<Texture> blueConstant = std::shared_ptr<Texture>(new ConstantTexture(Vec3(.25f, .25f, .75f)));
		std::shared_ptr<Texture> whiteConstant = std::shared_ptr<Texture>(new ConstantTexture(Vec3(.75f, .75f, .75f)));
		std::shared_ptr<Texture> fullWhiteConstant = std::shared_ptr<Texture>(new ConstantTexture(Vec3(1, 1, 1)));


		//Left
		std::shared_ptr<Shape> leftWallShape = std::shared_ptr<Shape>(new Sphere(1e5f, Vec3(1e5f + 1, 40.8f, 81.6f)));
		std::shared_ptr<Material> leftWallMaterial = std::shared_ptr<Material>(new DiffuseMaterial(redConstant));
		std::shared_ptr<Primitive> leftWall = std::shared_ptr<Primitive>(new GeometryPrimitive(leftWallShape, leftWallMaterial));
		scene->AddPrimitive(leftWall);

		//Right
		std::shared_ptr<Shape> rightWallShape = std::shared_ptr<Shape>(new Sphere(1e5f, Vec3(-1e5f + 99, 40.8f, 81.6f)));
		std::shared_ptr<Material> rightWallMaterial = std::shared_ptr<Material>(new DiffuseMaterial(blueConstant));
		std::shared_ptr<Primitive> rightWall = std::shared_ptr<Primitive>(new GeometryPrimitive(rightWallShape, rightWallMaterial));
		scene->AddPrimitive(rightWall);

		//Back
		std::shared_ptr<Shape> backWallShape = std::shared_ptr<Shape>(new Sphere(1e5f, Vec3(50, 40.8f, 1e5)));
		std::shared_ptr<Material> backWallMaterial = std::shared_ptr<Material>(new DiffuseMaterial(whiteConstant));
		std::shared_ptr<Primitive> backWall = std::shared_ptr<Primitive>(new GeometryPrimitive(backWallShape, backWallMaterial));
		scene->AddPrimitive(backWall);

		//Botom
		std::shared_ptr<Shape> botomShape = std::shared_ptr<Shape>(new Sphere(1e5f, Vec3(50, 1e5, 81.6f)));
		std::shared_ptr<Material> botomMaterial = std::shared_ptr<Material>(new DiffuseMaterial(whiteConstant));
		std::shared_ptr<Primitive> botom = std::shared_ptr<Primitive>(new GeometryPrimitive(botomShape, botomMaterial));
		scene->AddPrimitive(botom);

		//Top
		std::shared_ptr<Shape> topShape = std::shared_ptr<Shape>(new Sphere(1e5f, Vec3(50, -1e5f + 81.6f, 81.6f)));
		std::shared_ptr<Material> topMaterial = std::shared_ptr<Material>(new DiffuseMaterial(whiteConstant));
		std::shared_ptr<Primitive> top = std::shared_ptr<Primitive>(new GeometryPrimitive(topShape, topMaterial));
		scene->AddPrimitive(top);

		//Diffuse Ball1
		std::shared_ptr<Shape> diffuseBallShape1 = std::shared_ptr<Shape>(new Sphere(9.5f, Vec3(23, 0.0f, 98)));
		std::shared_ptr<Material> diffuseBallMaterial = std::shared_ptr<Material>(new DiffuseMaterial(fullWhiteConstant));
		std::shared_ptr<Primitive> diffuseBall1 = std::shared_ptr<Primitive>(new GeometryPrimitive(diffuseBallShape1, diffuseBallMaterial));
		scene->AddPrimitive(diffuseBall1);

		//Glass Ball1
		std::shared_ptr<Shape> glassBallShape1 = std::shared_ptr<Shape>(new Sphere(16.5f, Vec3(73, 26.5f, 78)));
		std::shared_ptr<Material> glassBallMaterial = std::shared_ptr<Material>(new GlassMaterial(fullWhiteConstant, fullWhiteConstant));
		std::shared_ptr<Primitive> glassBall1 = std::shared_ptr<Primitive>(new GeometryPrimitive(glassBallShape1, glassBallMaterial));
		scene->AddPrimitive(glassBall1);

		//Diffuse Ball2 new Sphere(7.0, Vec3(27, 16.5, 47), Vec3(), Vec3(.25, .25, .75)
		//std::shared_ptr<Shape> diffuseBallShape2 = std::shared_ptr<Shape>(new Sphere(7.0, Vec3(27, 16.5, 47)));
		//std::shared_ptr<Material> diffuseBallMaterial2 = std::shared_ptr<Material>(new DiffuseMaterial(blueConstant));
		//std::shared_ptr<Primitive> diffuseBall2 = std::shared_ptr<Primitive>(new GeometryPrimitive(diffuseBallShape2, diffuseBallMaterial2));
		//scene->AddPrimitive(diffuseBall2);


		//Glass Ball2 new Sphere(16.5f, Vec3(27, 16.5f, 47)
		std::shared_ptr<Shape> glassBallShape2 = std::shared_ptr<Shape>(new Sphere(16.5f, Vec3(27, 16.5f, 47)));
		std::shared_ptr<Primitive> glassBall2 = std::shared_ptr<Primitive>(new GeometryPrimitive(glassBallShape2, glassBallMaterial));
		scene->AddPrimitive(glassBall2);

		//Glass Ball3 new Sphere(9.5f, Vec3(53, 9.5f, 88)
		std::shared_ptr<Shape> glassBallShape3 = std::shared_ptr<Shape>(new Sphere(9.5f, Vec3(53, 9.5f, 88)));
		std::shared_ptr<Primitive> glassBall3 = std::shared_ptr<Primitive>(new GeometryPrimitive(glassBallShape3, glassBallMaterial));
		scene->AddPrimitive(glassBall3);

		//Light
		std::shared_ptr<Texture> lightTexture = std::shared_ptr<Texture>(new ConstantTexture(Vec3(0, 0, 0)));
		std::shared_ptr<Shape> lightShape = std::shared_ptr<Shape>(new Sphere(8.f, Vec3(50, 81.6f - 16.5f, 81.6f)));//Lite
		std::shared_ptr<Light> light0 = std::shared_ptr<Light>(new AreaLight(lightShape, Vec3(0.3f, 0.3f, 0.3f) * 100));
		std::shared_ptr<Material> lightMaterial = std::shared_ptr<Material>(new DiffuseMaterial(lightTexture));
		std::shared_ptr<Primitive> lightPrimitive = std::shared_ptr<Primitive>(new GeometryPrimitive(lightShape, lightMaterial, light0));
		scene->AddLight(lightPrimitive);

	}
};


class CornellBoxTriangle {
public:
	static void SetScene(std::shared_ptr<Scene>& scene) {
		//texture
		std::shared_ptr<Texture> redConstant = std::shared_ptr<Texture>(new ConstantTexture(Vec3(.75f, .25f, .25f)));
		std::shared_ptr<Texture> blueConstant = std::shared_ptr<Texture>(new ConstantTexture(Vec3(.25f, .25f, .75f)));
		std::shared_ptr<Texture> whiteConstant = std::shared_ptr<Texture>(new ConstantTexture(Vec3(.75f, .75f, .75f)));
		std::shared_ptr<Texture> fullWhiteConstant = std::shared_ptr<Texture>(new ConstantTexture(Vec3(1, 1, 1)));


		//Left
		std::shared_ptr<Shape> leftWallShape = std::shared_ptr<Shape>(new Sphere(1e5f, Vec3(1e5f + 1, 40.8f, 81.6f)));
		std::shared_ptr<Material> leftWallMaterial = std::shared_ptr<Material>(new DiffuseMaterial(redConstant));
		std::shared_ptr<Primitive> leftWall = std::shared_ptr<Primitive>(new GeometryPrimitive(leftWallShape, leftWallMaterial));
		scene->AddPrimitive(leftWall);

		//Right
		std::shared_ptr<Shape> rightWallShape = std::shared_ptr<Shape>(new Sphere(1e5f, Vec3(-1e5f + 99, 40.8f, 81.6f)));
		std::shared_ptr<Material> rightWallMaterial = std::shared_ptr<Material>(new DiffuseMaterial(blueConstant));
		std::shared_ptr<Primitive> rightWall = std::shared_ptr<Primitive>(new GeometryPrimitive(rightWallShape, rightWallMaterial));
		scene->AddPrimitive(rightWall);

		//Back
		std::shared_ptr<Shape> backWallShape = std::shared_ptr<Shape>(new Sphere(1e5f, Vec3(50, 40.8f, 1e5)));
		std::shared_ptr<Material> backWallMaterial = std::shared_ptr<Material>(new DiffuseMaterial(whiteConstant));
		std::shared_ptr<Primitive> backWall = std::shared_ptr<Primitive>(new GeometryPrimitive(backWallShape, backWallMaterial));
		scene->AddPrimitive(backWall);

		//Botom
		std::shared_ptr<Shape> botomShape = std::shared_ptr<Shape>(new Sphere(1e5f, Vec3(50, 1e5, 81.6f)));
		std::shared_ptr<Material> botomMaterial = std::shared_ptr<Material>(new DiffuseMaterial(whiteConstant));
		std::shared_ptr<Primitive> botom = std::shared_ptr<Primitive>(new GeometryPrimitive(botomShape, botomMaterial));
		scene->AddPrimitive(botom);

		//Top
		std::shared_ptr<Shape> topShape = std::shared_ptr<Shape>(new Sphere(1e5f, Vec3(50, -1e5f + 81.6f, 81.6f)));
		std::shared_ptr<Material> topMaterial = std::shared_ptr<Material>(new DiffuseMaterial(whiteConstant));
		std::shared_ptr<Primitive> top = std::shared_ptr<Primitive>(new GeometryPrimitive(topShape, topMaterial));
		scene->AddPrimitive(top);

		//Triangle face
		//Transform transformTriangle = Transform::Translate(Vec3(50, 52, 70.6f));
		//Vec3 p0 = transformTriangle(Vec3(0, 5, -5));
		//Vec3 p1 = transformTriangle(Vec3(0, 5, 5));
		//Vec3 p2 = transformTriangle(Vec3(0, -5, 5));
		Transform transformTriangle = Transform::Translate(Vec3(27, 16.5, 47)) * Transform::RotateY(30);
		Vec3 p0 = transformTriangle(Vec3(10, 20, 0));
		Vec3 p1 = transformTriangle(Vec3(-10, 0, 0));
		Vec3 p2 = transformTriangle(Vec3(10, 0, 0));
		Vec3 faceNormal = Cross(p1 - p0, p2 - p0).Norm();
		std::shared_ptr<Shape> triangleShape =
			std::shared_ptr<Shape>(new Triangle(p0, p1, p2, faceNormal));
		std::shared_ptr<Material> triangleMaterial = std::shared_ptr<Material>(new DiffuseMaterial(blueConstant));
		std::shared_ptr<Primitive> triangle = std::shared_ptr<Primitive>(new GeometryPrimitive(triangleShape, triangleMaterial));
		scene->AddPrimitive(triangle);

		////Diffuse Ball1
		//std::shared_ptr<Shape> diffuseBallShape1 = std::shared_ptr<Shape>(new Sphere(9.5f, Vec3(23, 0.0f, 98)));
		//std::shared_ptr<Material> diffuseBallMaterial = std::shared_ptr<Material>(new DiffuseMaterial(fullWhiteConstant));
		//std::shared_ptr<Primitive> diffuseBall1 = std::shared_ptr<Primitive>(new GeometryPrimitive(diffuseBallShape1, diffuseBallMaterial));
		//scene->AddPrimitive(diffuseBall1);

		////Glass Ball1
		//std::shared_ptr<Shape> glassBallShape1 = std::shared_ptr<Shape>(new Sphere(16.5f, Vec3(73, 26.5f, 78)));
		//std::shared_ptr<Material> glassBallMaterial = std::shared_ptr<Material>(new GlassMaterial(fullWhiteConstant, fullWhiteConstant));
		//std::shared_ptr<Primitive> glassBall1 = std::shared_ptr<Primitive>(new GeometryPrimitive(glassBallShape1, glassBallMaterial));
		//scene->AddPrimitive(glassBall1);

		////Diffuse Ball2 new Sphere(7.0, Vec3(27, 16.5, 47), Vec3(), Vec3(.25, .25, .75)
		////std::shared_ptr<Shape> diffuseBallShape2 = std::shared_ptr<Shape>(new Sphere(7.0, Vec3(27, 16.5, 47)));
		////std::shared_ptr<Material> diffuseBallMaterial2 = std::shared_ptr<Material>(new DiffuseMaterial(blueConstant));
		////std::shared_ptr<Primitive> diffuseBall2 = std::shared_ptr<Primitive>(new GeometryPrimitive(diffuseBallShape2, diffuseBallMaterial2));
		////scene->AddPrimitive(diffuseBall2);


		////Glass Ball2 new Sphere(16.5f, Vec3(27, 16.5f, 47)
		//std::shared_ptr<Shape> glassBallShape2 = std::shared_ptr<Shape>(new Sphere(16.5f, Vec3(27, 16.5f, 47)));
		//std::shared_ptr<Primitive> glassBall2 = std::shared_ptr<Primitive>(new GeometryPrimitive(glassBallShape2, glassBallMaterial));
		//scene->AddPrimitive(glassBall2);

		////Glass Ball3 new Sphere(9.5f, Vec3(53, 9.5f, 88)
		//std::shared_ptr<Shape> glassBallShape3 = std::shared_ptr<Shape>(new Sphere(9.5f, Vec3(53, 9.5f, 88)));
		//std::shared_ptr<Primitive> glassBall3 = std::shared_ptr<Primitive>(new GeometryPrimitive(glassBallShape3, glassBallMaterial));
		//scene->AddPrimitive(glassBall3);

		//Light
		//std::shared_ptr<Texture> lightTexture = std::shared_ptr<Texture>(new ConstantTexture(Vec3(0, 0, 0)));
		//std::shared_ptr<Shape> lightShape = std::shared_ptr<Shape>(new Sphere(8.f, Vec3(50, 81.6f - 16.5f, 81.6f)));//Lite
		//std::shared_ptr<Light> light0 = std::shared_ptr<Light>(new AreaLight(lightShape, Vec3(0.3f, 0.3f, 0.3f) * 100));
		//std::shared_ptr<Material> lightMaterial = std::shared_ptr<Material>(new DiffuseMaterial(lightTexture));
		//std::shared_ptr<Primitive> lightPrimitive = std::shared_ptr<Primitive>(new GeometryPrimitive(lightShape, lightMaterial, light0));
		//scene->AddLight(lightPrimitive);

		//Light
		Transform lightTransform = Transform::Translate(Vec3(50, 81.6f - 16.5f, 81.6f));
		Vec3 lightP0 = lightTransform(Vec3(6, 0, 6));
		Vec3 lightP1 = lightTransform(Vec3(6, 0, -6));
		Vec3 lightP2 = lightTransform(Vec3(-6, 0, 6));
		Vec3 lightP3 = lightTransform(Vec3(-6, 0, -6));
		Vec3 lightNormal = Vec3(0, -1, 0);
		std::shared_ptr<Texture> lightTexture = std::shared_ptr<Texture>(new ConstantTexture(Vec3(0, 0, 0)));
		std::shared_ptr<Shape> lightTriangleShape0 = 
			std::shared_ptr<Shape>(new Triangle(lightP0, lightP1, lightP2, lightNormal));
		std::shared_ptr<Shape> lightTriangleShape1 =
			std::shared_ptr<Shape>(new Triangle(lightP1, lightP3, lightP2, lightNormal));
		std::shared_ptr<Light> light0 = std::shared_ptr<Light>(new AreaLight(lightTriangleShape0, Vec3(0.3f, 0.3f, 0.3f) * 100));
		std::shared_ptr<Light> light1 = std::shared_ptr<Light>(new AreaLight(lightTriangleShape1, Vec3(0.3f, 0.3f, 0.3f) * 100));
		std::shared_ptr<Material> lightMaterial = std::shared_ptr<Material>(new DiffuseMaterial(lightTexture));
		std::shared_ptr<Primitive> lightPrimitive0 = 
			std::shared_ptr<Primitive>(new GeometryPrimitive(lightTriangleShape0, lightMaterial, light0));
		std::shared_ptr<Primitive> lightPrimitive1 = 
			std::shared_ptr<Primitive>(new GeometryPrimitive(lightTriangleShape1, lightMaterial, light1));
		scene->AddLight(lightPrimitive0);
		//scene->AddLight(lightPrimitive1);
	}
};
NAMESPACE_END