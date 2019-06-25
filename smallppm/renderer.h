#pragma once

#include "integrator.h"
#include "scene.h"

class Renderer {
public:
	Renderer(const std::shared_ptr<Scene> &pScene, const std::shared_ptr<Integrator> &pIntegrator,
		const std::shared_ptr<Film> &pFilm) : scene(pScene), integrator(pIntegrator), film(pFilm) {

	}
	void Render() {
		integrator->Render(*scene);
		film->SaveImage();
	}
private:
	std::shared_ptr<Scene> scene;
	std::shared_ptr<Integrator> integrator;
	std::shared_ptr<Film> film;
};
