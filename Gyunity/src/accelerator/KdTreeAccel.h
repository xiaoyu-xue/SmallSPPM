#pragma once

#include "common/Core.h"
#include "Accelerator.h"

GY_NAMESPACE_BEGIN

struct KdAccelNode;
struct BoundEdge;
class KdTreeAccel : public Accelerator{
public:
	// KdTreeAccel Public Methods
	KdTreeAccel(std::vector<std::shared_ptr<Primitive>> p,
		int isectCost = 80, int traversalCost = 1,
		real emptyBonus = 0.5, int maxPrims = 1, int maxDepth = -1);
	AABB WorldBound() const { return bounds; }
	~KdTreeAccel();
	virtual bool Intersect(const Ray& r, Intersection* isect) const override;
	virtual bool Intersect(const Ray& r) const override;
	//Todo
	void SetPrimitives(const std::vector<std::shared_ptr<Primitive>>& a) override {}
private:
	// KdTreeAccel Private Methods
	void buildTree(int nodeNum, const AABB& bounds,
		const std::vector<AABB>& primBounds, int* primNums,
		int nprims, int depth,
		const std::unique_ptr<BoundEdge[]> edges[3], int* prims0,
		int* prims1, int badRefines = 0);

	// KdTreeAccel Private Data
	const int isectCost, traversalCost, maxPrims;
	const real emptyBonus;
	std::vector<std::shared_ptr<Primitive>> primitives;
	std::vector<int> primitiveIndices;
	KdAccelNode* nodes;
	int nAllocedNodes, nextFreeNode;
	AABB bounds;
};

struct KdToDo {
	const KdAccelNode* node;
	real tMin, tMax;
};

GY_NAMESPACE_END