#pragma once

#include "utils.h"
#include "linagl.h"
#include <vector>
#include "AABB.h"
#include "hitpoint.h"
#include <mutex>
#include "threading.h"

NAMESPACE_BEGIN

template<typename T>
class HashGrid {
public:
	HashGrid(real irad = 0.8) : searchRadius(irad) {}

	void Initialize(real irad) {
		searchRadius = irad;
		hpbbox.Reset();
	}

	inline void SetSearchRadius(real radius) {
		searchRadius = radius;
	}

	// spatial hash function
	inline uint32 hash(const int ix, const int iy, const int iz) {
		return (uint32)((ix * 73856093) ^ (iy * 19349663) ^ (iz * 83492791)) % hashNum;
	}

	inline void AddPoint(const std::pair<Vec, T> &point) {
		hitPoints.push_back(point);
		hpbbox.Fit(point.first);
	}

	inline void AddPoint(std::pair<Vec, T> &&point) {
		hitPoints.push_back(point);
		hpbbox.Fit(point.first);
	}

	void BuildHashGrid(real radius = 0.8) {
		// find the bounding box of all the measurement points
		//hpbbox.Reset();
		//for (const std::pair<Vec, T> &hp : hitPoints) {
		//	hpbbox.Fit(hp.first);
		//}

		real searchRadius = radius;
		// heuristic for initial radius
		Vec ssize = hpbbox.maxPoint - hpbbox.minPoint;
		// determine hash table size
		// we now find the bounding box of all the measurement points inflated by the initial radius
		hpbbox.Reset();
		int64 vphoton = 0;

		for (const std::pair<Vec, T> &hp : hitPoints) {
			vphoton++;
			hpbbox.Fit(hp.first - searchRadius);
			hpbbox.Fit(hp.first + searchRadius);
		}

		// make each grid cell two times larger than the initial radius
		invHashCellSize = 1.0 / (searchRadius * 2.0);
		hashNum = vphoton;

		// build the hash table

		hashGrid.resize(hashNum);
		hashGridSpinlocks.resize(hashNum);


		int hitPointNum = (int)hitPoints.size();
		ParallelFor(0, hitPointNum, [&](int i) {
			std::pair<Vec, T> &hp = hitPoints[i];
				Vec BMin = ((hp.first - searchRadius) - hpbbox.minPoint) * invHashCellSize;
				Vec BMax = ((hp.first + searchRadius) - hpbbox.minPoint) * invHashCellSize;
				for (int iz = abs(int(BMin.z)); iz <= abs(int(BMax.z)); iz++)
				{
					for (int iy = abs(int(BMin.y)); iy <= abs(int(BMax.y)); iy++)
					{
						for (int ix = abs(int(BMin.x)); ix <= abs(int(BMax.x)); ix++)
						{
							int hv = hash(ix, iy, iz);
							std::lock_guard<Spinlock> lock(hashGridSpinlocks[hv]);
							hashGrid[hv].push_back(hp.second);
						}
					}
				}
		});
	}

	std::vector<T>& GetGrid(const Vec &point) {
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
		//std::vector<std::vector<T>>().swap(hashGrid);
		//std::vector<std::pair<Vec, T>>().swap(hitPoints);
		hashGrid.clear();
		hitPoints.clear();
		hpbbox.Reset();
	}

private:
	std::vector<std::vector<T>> hashGrid;
	std::vector<std::pair<Vec, T>> hitPoints;
	std::vector<Spinlock> hashGridSpinlocks;
	int64 hashNum, pixelIndex;
	real invHashCellSize;
	real searchRadius;
	AABB hpbbox;
};


NAMESPACE_END