#pragma once

#include "utils.h"

NAMESPACE_BEGIN

class Filter
{
protected:
	real radius;

public:
	Filter(const real rad)
		: radius(rad)
	{
	}
	virtual ~Filter() {}

	const real GetRadius() const
	{
		return radius;
	}
	virtual real Evaluate(const real dx, const real dy) const = 0;
};

class BoxFilter : public Filter
{
public:
	BoxFilter()
		: Filter(0.5)
	{
	}

public:
	real Evaluate(const real dx, const real dy) const {
		return 1.0;
	}
};

NAMESPACE_END