#pragma once

#include "utils.h"
#include "linagl.h"
#include <vector>
//#include "threading.h"
#include "AABB.h"
#include "hitpoint.h"
#include <mutex>
#include "threading.h"

class HashGrid {
public:
	HashGrid(real initialRadius = 1.0) {
		irad = initialRadius;
	}

	void Initialize(real initialRadius) {
		irad = initialRadius;
	}

	// spatial hash function
	inline uint32 hash(const int ix, const int iy, const int iz) {
		return (uint32)((ix * 73856093) ^ (iy * 19349663) ^ (iz * 83492791)) % hashNum;
	}

	void BuildHashGrid(std::vector<HPoint> &hitPoints) {
		// find the bounding box of all the measurement points
		hpbbox.reset();
		for (const HPoint &hp : hitPoints) {
			if (hp.used) hpbbox.fit(hp.pos);
		}

		// heuristic for initial radius
		Vec ssize = hpbbox.maxPoint - hpbbox.minPoint;
		//real irad = ((ssize.x + ssize.y + ssize.z) / 3.0) / ((w + h) / 2.0) * 2.0;
		real irad = 0.8;
		// determine hash table size
		// we now find the bounding box of all the measurement points inflated by the initial radius
		hpbbox.reset();
		int vphoton = 0;

		for (const HPoint &hp : hitPoints) {
			if (!hp.used) continue;
			vphoton++;
			hpbbox.fit(hp.pos - irad);
			hpbbox.fit(hp.pos + irad);
		}

		// make each grid cell two times larger than the initial radius
		invHashCellSize = 1.0 / (irad * 2.0);
		hashNum = vphoton;

		// build the hash table

		hashGrid.resize(hashNum);
		hashGridSpinlocks.resize(hashNum);


		//for (HPoint &hp : hitPoints) {
		int hitPointNum = (int)hitPoints.size();
//#pragma omp parallel for schedule(guided)
		//for (int i = 0; i < hitPointNum; ++i) {
		ParallelFor(0, hitPointNum, [&](int i) {
			HPoint &hp = hitPoints[i];
			if (hp.used) {
				Vec BMin = ((hp.pos - irad) - hpbbox.minPoint) * invHashCellSize;
				Vec BMax = ((hp.pos + irad) - hpbbox.minPoint) * invHashCellSize;
				for (int iz = abs(int(BMin.z)); iz <= abs(int(BMax.z)); iz++)
				{
					for (int iy = abs(int(BMin.y)); iy <= abs(int(BMax.y)); iy++)
					{
						for (int ix = abs(int(BMin.x)); ix <= abs(int(BMax.x)); ix++)
						{
							int hv = hash(ix, iy, iz);
							std::lock_guard<Spinlock> lock(hashGridSpinlocks[hv]);
							hashGrid[hv].push_back(&hp);
						}
					}
				}
			}
		});
		//}
	}

	std::vector<HPoint*>& GetGrid(const Vec &point) {
		// photon ray
		// find neighboring measurement points and accumulate flux via progressive density estimation
		Vec hh = (point - hpbbox.minPoint) * invHashCellSize;
		int ix = std::abs(int(hh.x)), iy = std::abs(int(hh.y)), iz = std::abs(int(hh.z));
		// strictly speaking, we should use #pragma omp critical here.
		// it usually works without an artifact due to the fact that photons are 
		// rarely accumulated to the same measurement points at the same time (especially with QMC).
		// it is also significantly faster.
		return hashGrid[hash(ix, iy, iz)];
	}


	void ClearHashGrid() {
		for (auto &e : hashGrid) {
			std::vector<HPoint*>().swap(e);
		}
	}

private:
	std::vector<std::vector<HPoint*>> hashGrid;
	std::vector<Spinlock> hashGridSpinlocks;
	int64 hashNum, pixelIndex;
	real invHashCellSize;
	real irad;
	AABB hpbbox;
};