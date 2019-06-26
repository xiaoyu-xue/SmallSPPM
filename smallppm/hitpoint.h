#pragma once

#include "utils.h"
#include "linagl.h"

NAMESPACE_BEGIN

struct HPoint {
	Vec importance, pos, nrm, flux, outDir;
	real r2;
	int64 n; // n = N / ALPHA in the paper
	int64 pix;
	bool used;
};

NAMESPACE_END