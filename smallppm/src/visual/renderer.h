#pragma once

#include "integrator/Integrator.h"
#include "Scene.h"

NAMESPACE_BEGIN

class Renderer {
public:
	Renderer(const std::shared_ptr<Scene> &pScene, const std::shared_ptr<Camera> &pCamera, 
		const std::shared_ptr<Integrator> &pIntegrator, const std::shared_ptr<Film> &pFilm) :
		scene(pScene), integrator(pIntegrator), film(pFilm), camera(pCamera) {

	}
	void Render() {
		integrator->Render(*scene, *camera);
		film->SaveImage();
	}
private:
	std::shared_ptr<Scene> scene;
	std::shared_ptr<Integrator> integrator;
	std::shared_ptr<Film> film;
	std::shared_ptr<Camera> camera;
	
};

NAMESPACE_END