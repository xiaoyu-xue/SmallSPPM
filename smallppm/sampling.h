#pragma once
#include "utils.h"
#include <string.h>
#include <vector>
#include <algorithm>
#include "scalar.h"
#include "def.h"
#include "linagl.h"

NAMESPACE_BEGIN

struct Distribution1D
{
public:
	Distribution1D(const real *f, int n)
	{
		count = n;
		func = new real[n];
		memcpy(func, f, n * sizeof(real));
		cdf = new real[n + 1];

		// Compute integral of step function at $x_i$
		cdf[0] = 0.;
		for (int i = 1; i < count + 1; ++i)
			cdf[i] = cdf[i - 1] + func[i - 1] / n;

		// Transform step function integral into CDF
		funcInt = cdf[count];
		if (funcInt == 0.f)
		{
			for (int i = 1; i < n + 1; ++i)
				cdf[i] = real(i) / real(n);
		}
		else
		{
			for (int i = 1; i < n + 1; ++i)
				cdf[i] /= funcInt;
		}
	}

	~Distribution1D()
	{
		delete[] func;
		delete[] cdf;
	}

	real SampleContinuous(real u, real *pdf, int *off = NULL) const
	{
		// Find surrounding CDF segments and _offset_
		real *ptr = std::upper_bound(cdf, cdf + count + 1, u);
		int offset = Clamp(int(ptr - cdf - 1), 0, count - 1);
		if (off) *off = offset;
		ASSERT(offset < count);
		ASSERT(u >= cdf[offset] && (u < cdf[offset + 1] || u == 1));

		// Fix the case when func ends with zeros
		if (cdf[offset] == cdf[offset + 1])
		{
			ASSERT(u == 1.0f);

			do { offset--; } while (cdf[offset] == cdf[offset + 1] && offset > 0);

			ASSERT(cdf[offset] != cdf[offset + 1]);
		}

		// Compute offset along CDF segment
		real du = (u - cdf[offset]) / (cdf[offset + 1] - cdf[offset]);
		ASSERT(!std::isnan(du));

		// Compute PDF for sampled offset
		if (pdf) *pdf = func[offset] / funcInt;
		ASSERT(func[offset] > 0);

		// Return $x\in{}[0,1]$ corresponding to sample
		return (offset + du) / count;
	}

	int SampleDiscrete(real u, real *pdf) const
	{
		// Find surrounding CDF segments and _offset_
		real *ptr = std::upper_bound(cdf, cdf + count + 1, u);
		int offset = std::max(0, int(ptr - cdf - 1));
		ASSERT(offset < count);
		ASSERT(u >= cdf[offset] && u < cdf[offset + 1]);
		if (pdf) *pdf = func[offset] / (funcInt * count);
		return offset;
	}

private:
	friend struct Distribution2D;
	real *func, *cdf;
	real funcInt;
	int count;
};

struct Distribution2D
{
public:
	Distribution2D(const real *func, int nu, int nv)
	{
		pConditionalV.reserve(nv);
		for (int v = 0; v < nv; ++v)
		{
			// Compute conditional sampling distribution for $\tilde{v}$
			pConditionalV.push_back(new Distribution1D(&func[v*nu], nu));
		}

		// Compute marginal sampling distribution $p[\tilde{v}]$
		std::vector<real> marginalFunc;
		marginalFunc.reserve(nv);
		for (int v = 0; v < nv; ++v)
			marginalFunc.push_back(pConditionalV[v]->funcInt);
		pMarginal = new Distribution1D(&marginalFunc[0], nv);
	}

	~Distribution2D()
	{
		delete pMarginal;
		for (uint32_t i = 0; i < pConditionalV.size(); ++i)
			delete pConditionalV[i];
	}

	void SampleContinuous(real u0, real u1, real uv[2], real *pdf) const
	{
		real pdfs[2];
		int v;
		uv[1] = pMarginal->SampleContinuous(u1, &pdfs[1], &v);
		uv[0] = pConditionalV[v]->SampleContinuous(u0, &pdfs[0]);
		*pdf = pdfs[0] * pdfs[1];
	}

	real Pdf(real u, real v) const
	{
		int iu = Clamp((int)(u * pConditionalV[0]->count), 0,
			pConditionalV[0]->count - 1);
		int iv = Clamp((int)(v * pMarginal->count), 0,
			pMarginal->count - 1);
		if (pConditionalV[iv]->funcInt * pMarginal->funcInt == 0.f) return 0.f;
		return (pConditionalV[iv]->func[iu] * pMarginal->func[iv]) /
			(pConditionalV[iv]->funcInt * pMarginal->funcInt);
	}

private:
	std::vector<Distribution1D *> pConditionalV;
	Distribution1D *pMarginal;
};



Vec2 ConcentricSampleDisk(const Vec2 &u);

Vec3 UniformSampleSphere(const Vec2 &u);

Vec3 CosineSampleHemisphere(const Vec2 &u);

real UniformSpherePdf();

Vec2 UniformSampleTriangle(const Vec2 &u);

real BalanceHeuristic(int nf, real fPdf, int ng, real gPdf);

real CosineHemispherePdf(real cosTheta);


NAMESPACE_END