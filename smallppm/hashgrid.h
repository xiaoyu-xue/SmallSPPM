#pragma once

#include "utils.h"
#include "linagl.h"
#include <vector>
#include "AABB.h"
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


	// spatial hash function
	inline uint32 hash(const int ix, const int iy, const int iz) {
		return (uint32)((ix * 73856093) ^ (iy * 19349663) ^ (iz * 83492791)) % hashNum;
	}

	inline void AddPoint(const std::pair<Vec3, T> &point) {
		hitPoints.push_back(point);
		hpbbox.Fit(point.first);
	}

	inline void AddPoint(std::pair<Vec3, T> &&point) {
		hitPoints.push_back(point);
		hpbbox.Fit(point.first);
	}

	inline void AddPoint(const std::pair<Vec3, T> &point, real radius) {
		hitPoints.push_back(point);
		radii.push_back(radius);
		hpbbox.Fit(point.first);
	}

	inline void AddPoint(std::pair<Vec3, T> &&point, real radius) {
		hitPoints.push_back(point);
		radii.push_back(radius);
		hpbbox.Fit(point.first);
	}

	void BuildHashGrid(real maxRadius = 0.8f) {
		// find the bounding box of all the measurement points
		//hpbbox.Reset();
		//for (const std::pair<Vec3, T> &hp : hitPoints) {
		//	hpbbox.Fit(hp.first);
		//}

		real searchRadius = maxRadius;
		// heuristic for initial radius
		//Vec3 ssize = hpbbox.maxPoint - hpbbox.minPoint;
		// determine hash table size
		// we now find the bounding box of all the measurement points inflated by the initial radius
		hpbbox.Reset();
		int64 vphoton = 0;

		for (const std::pair<Vec3, T> &hp : hitPoints) {
			vphoton++;
			hpbbox.Fit(hp.first - searchRadius);
			hpbbox.Fit(hp.first + searchRadius);
		}

		// make each grid cell two times larger than the initial radius
		invHashCellSize = 1.f / (searchRadius * 2.f);
		hashNum = vphoton;

		// build the hash table

		hashGrid.resize(hashNum);
		hashGridSpinlocks.resize(hashNum);

		int hitPointNum = (int)hitPoints.size();
		ParallelFor(0, hitPointNum, [&](int i) {
			std::pair<Vec3, T> &hp = hitPoints[i];
				Vec3 BMin = ((hp.first - radii[i]) - hpbbox.minPoint) * invHashCellSize;
				Vec3 BMax = ((hp.first + radii[i]) - hpbbox.minPoint) * invHashCellSize;
				for (int iz = std::abs(int(BMin.z)); iz <= std::abs(int(BMax.z)); iz++)
				{
					for (int iy = std::abs(int(BMin.y)); iy <= std::abs(int(BMax.y)); iy++)
					{
						for (int ix = std::abs(int(BMin.x)); ix <= std::abs(int(BMax.x)); ix++)
						{
							int hv = hash(ix, iy, iz);
							std::lock_guard<Spinlock> lock(hashGridSpinlocks[hv]);
							hashGrid[hv].push_back(hp.second);
						}
					}
				}
		});
	}

	std::vector<T>& GetGrid(const Vec3 &point) {
		Vec3 hh = (point - hpbbox.minPoint) * invHashCellSize;
		int ix = std::abs(int(hh.x)), iy = std::abs(int(hh.y)), iz = std::abs(int(hh.z));
		return hashGrid[hash(ix, iy, iz)];
	}


	void ClearHashGrid() {
		//std::vector<std::vector<T>>().swap(hashGrid);
		//std::vector<std::pair<Vec3, T>>().swap(hitPoints);
		hashGrid.clear();
		hitPoints.clear();
		radii.clear();
		hpbbox.Reset();
	}

private:
	std::vector<std::vector<T>> hashGrid;
	std::vector<std::pair<Vec3, T>> hitPoints;
	std::vector<real> radii;
	std::vector<Spinlock> hashGridSpinlocks;
	int64 hashNum;
	real invHashCellSize;
	real searchRadius;
	AABB hpbbox;
};


NAMESPACE_END