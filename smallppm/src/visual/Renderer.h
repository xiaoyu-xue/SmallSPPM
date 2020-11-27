#pragma once

#include "integrator/Integrator.h"
#include "Scene.h"

NAMESPACE_BEGIN

class Renderer 
{
private:
	std::shared_ptr<Scene> mpScene;
	std::shared_ptr<Integrator> mpIntegrator;
	std::shared_ptr<Film> mpFilm;
	std::shared_ptr<Camera> mpCamera;
public:
	Renderer(const std::shared_ptr<Scene> &pScene, const std::shared_ptr<Camera> &pCamera, 
		const std::shared_ptr<Integrator> &pIntegrator, const std::shared_ptr<Film> &pFilm) 
		: mpScene(pScene), mpIntegrator(pIntegrator), mpFilm(pFilm), mpCamera(pCamera) 
	{

	}
	void Render() 
	{
		mpIntegrator->Render(*mpScene, *mpCamera);
		mpFilm->SaveImage();
	}
};

NAMESPACE_END