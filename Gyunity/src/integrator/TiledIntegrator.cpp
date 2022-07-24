#include "TiledIntegrator.h"
#include "sampler/Sampler.h"
#include "sampler/SamplerEnum.h"

GYT_NAMESPACE_BEGIN

void TiledIntegrator::Render(const Scene& scene, const Camera& camera) {

	fprintf(stderr, "Rendering ...\n");

	std::shared_ptr<Film> film = camera.GetFilm();
	int resX = film->mResX, resY = film->mResY;
	const int tileSize = 32;
	std::vector<Tile> tiles;
	int tileIndex = 0;
	for (int y = 0; y < resY; y += tileSize) {
		for (int x = 0; x < resX; x += tileSize) {
			Tile tile;
			tile.minX = x;
			tile.minY = y;
			tile.maxX = std::min(resX, x + tileSize);
			tile.maxY = std::min(resY, y + tileSize);
			tiles.push_back(tile);
		}
	}
	std::atomic<int64> workDone = 0;
	const int tileNum = tiles.size();
	ParallelFor(0, tileNum, [&](int i) {
		MemoryPool arena;
		const Tile& tile = tiles[i];
		std::shared_ptr<Sampler> tileSampler = mpSampler->Clone(i);
		for (int y = tile.minY; y < tile.maxY; ++y) {
			for (int x = tile.minX; x < tile.maxX; ++x) {
				for (int s = 0; s < spp; ++s) {
					uint64 instance = mpSamplerEnum->GetIndex(s, x, y);
					RandomStateSequence rand(tileSampler, instance);
					real u = mpSamplerEnum->SampleX(x, rand());
					real v = mpSamplerEnum->SampleY(y, rand());
					Vec2 pixelSample(u, v);
					Ray ray = camera.GenerateRay(x, y, pixelSample);
					Vec3 L = Li(ray, scene, rand, arena);
					film->AddSample(x + u, y + v, L);
					arena.Reset();
				}
			}
		}
		workDone += 1;
		real percentage = 100.f * workDone / tileNum;
		fprintf(stderr, "\rPercentage: %5.2f%%", percentage);
	});
}

GYT_NAMESPACE_END