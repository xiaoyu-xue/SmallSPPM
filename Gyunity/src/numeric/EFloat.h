#pragma once

#include "common/Core.h"
#include "NumericUtils.h"
#include <iostream>

GY_NAMESPACE_BEGIN

class EReal {
public:
	// EReal Public Methods
	EReal() {}
	EReal(real v, real err = 0.0) : v(v) {
		if (err == 0.)
			low = high = v;
		else {
			// Compute conservative bounds by rounding the endpoints away
			// from the middle. Note that this will be over-conservative in
			// cases where v-err or v+err are exactly representable in
			// floating-point, but it's probably not worth the trouble of
			// checking this case.
			low = NextFloatDown(v - err);
			high = NextFloatUp(v + err);
		}
		// Store high precision reference value in _EReal_
#ifndef NDEBUG
		vPrecise = v;
		//Check();
#endif  // NDEBUG
		}
#ifndef NDEBUG
	EReal(real v, long double lD, real err) : EReal(v, err) {
		vPrecise = lD;
	}
#endif  // DEBUG

	EReal operator+(EReal ef) const {
		EReal r;
		r.v = v + ef.v;
#ifndef NDEBUG
		r.vPrecise = vPrecise + ef.vPrecise;
#endif  // DEBUG
		// Interval arithemetic addition, with the result rounded away from
		// the value r.v in order to be conservative.
		r.low = NextFloatDown(LowerBound() + ef.LowerBound());
		r.high = NextFloatUp(UpperBound() + ef.UpperBound());
		//r.Check();
		return r;
	}

	explicit operator real() const { return v; }
	real GetAbsoluteError() const { return high - low; }
	real UpperBound() const { return high; }
	real LowerBound() const { return low; }

#ifndef NDEBUG
	real GetRelativeError() const {
		return std::abs((vPrecise - v) / vPrecise);
	}
	long double PreciseValue() const { return vPrecise; }
#endif

	EReal operator-(EReal ef) const {
		EReal r;
		r.v = v - ef.v;
#ifndef NDEBUG
		r.vPrecise = vPrecise - ef.vPrecise;
#endif
		r.low = NextFloatDown(LowerBound() - ef.UpperBound());
		r.high = NextFloatUp(UpperBound() - ef.LowerBound());
		//r.Check();
		return r;
	}
	EReal operator*(EReal ef) const {
		EReal r;
		r.v = v * ef.v;
#ifndef NDEBUG
		r.vPrecise = vPrecise * ef.vPrecise;
#endif
		real prod[4] = {
			LowerBound() * ef.LowerBound(), UpperBound() * ef.LowerBound(),
			LowerBound() * ef.UpperBound(), UpperBound() * ef.UpperBound() };
		r.low = NextFloatDown(
			std::min(std::min(prod[0], prod[1]), std::min(prod[2], prod[3])));
		r.high = NextFloatUp(
			std::max(std::max(prod[0], prod[1]), std::max(prod[2], prod[3])));
		//r.Check();
		return r;
	}

	EReal operator/(EReal ef) const {
		EReal r;
		r.v = v / ef.v;
#ifndef NDEBUG
		r.vPrecise = vPrecise / ef.vPrecise;
#endif
		if (ef.low < 0 && ef.high > 0) {
			// Bah. The interval we're dividing by straddles zero, so just
			// return an interval of everything.
			r.low = -Infinity;
			r.high = Infinity;
		}
		else {
			real div[4] = {
				LowerBound() / ef.LowerBound(), UpperBound() / ef.LowerBound(),
				LowerBound() / ef.UpperBound(), UpperBound() / ef.UpperBound() };
			r.low = NextFloatDown(
				std::min(std::min(div[0], div[1]), std::min(div[2], div[3])));
			r.high = NextFloatUp(
				std::max(std::max(div[0], div[1]), std::max(div[2], div[3])));
		}
		//r.Check();
		return r;
	}
	EReal operator-() const {
		EReal r;
		r.v = -v;
#ifndef NDEBUG
		r.vPrecise = -vPrecise;
#endif
		r.low = -high;
		r.high = -low;
		//r.Check();
		return r;
	}
	FORCE_INLINE bool operator==(EReal fe) const { return v == fe.v; }
//	FORCE_INLINE void Check() const {
//		if (!std::isinf(low) && !std::isnan(low) && !std::isinf(high) &&
//			!std::isnan(high))
//			CHECK_LE(low, high);
//#ifndef NDEBUG
//		if (!std::isinf(v) && !std::isnan(v)) {
//			CHECK_LE(LowerBound(), vPrecise);
//			CHECK_LE(vPrecise, UpperBound());
//		}
//#endif
//	}
	EReal(const EReal &ef) {
		//ef.Check();
		v = ef.v;
		low = ef.low;
		high = ef.high;
#ifndef NDEBUG
		vPrecise = ef.vPrecise;
#endif
	}
	EReal &operator=(const EReal &ef) {
		//ef.Check();
		if (&ef != this) {
			v = ef.v;
			low = ef.low;
			high = ef.high;
#ifndef NDEBUG
			vPrecise = ef.vPrecise;
#endif
		}
		return *this;
	}

	friend std::ostream &operator<<(std::ostream &os, const EReal &ef) {
		printf("v=%f (%a) - [%f, %f]",	ef.v, ef.v, ef.low, ef.high);
		//os << "v=" << ef.v << " - "
		//	<< "[" << ef.low << ", " << ef.high << "] ";
#ifndef NDEBUG
		printf(", precise=%.30Lf", ef.vPrecise);
		//os << ", precise=" << ef.vPrecise;
#endif // !NDEBUG
		return os;
	}

private:
	// EReal Private Data
	real v, low, high;
#ifndef NDEBUG
	long double vPrecise;
#endif  // NDEBUG
	friend FORCE_INLINE EReal Sqrt(EReal fe);
	friend FORCE_INLINE EReal Abs(EReal fe);
	friend FORCE_INLINE bool Quadratic(EReal A, EReal B, EReal C, EReal *t0, EReal *t1);
};

// EReal FORCE_INLINE Functions
FORCE_INLINE EReal operator*(float f, EReal fe) { return EReal(f) * fe; }

FORCE_INLINE EReal operator/(float f, EReal fe) { return EReal(f) / fe; }

FORCE_INLINE EReal operator+(float f, EReal fe) { return EReal(f) + fe; }

FORCE_INLINE EReal operator-(float f, EReal fe) { return EReal(f) - fe; }

FORCE_INLINE EReal Sqrt(EReal fe) {
	EReal r;
	r.v = std::sqrt(fe.v);
#ifndef NDEBUG
	r.vPrecise = std::sqrt(fe.vPrecise);
#endif
	r.low = NextFloatDown(std::sqrt(fe.low));
	r.high = NextFloatUp(std::sqrt(fe.high));
	//r.Check();
	return r;
}

FORCE_INLINE EReal Abs(EReal fe) {
	if (fe.low >= 0)
		// The entire interval is greater than zero, so we're all set.
		return fe;
	else if (fe.high <= 0) {
		// The entire interval is less than zero.
		EReal r;
		r.v = -fe.v;
#ifndef NDEBUG
		r.vPrecise = -fe.vPrecise;
#endif
		r.low = -fe.high;
		r.high = -fe.low;
		//r.Check();
		return r;
	}
	else {
		// The interval straddles zero.
		EReal r;
		r.v = std::abs(fe.v);
#ifndef NDEBUG
		r.vPrecise = std::abs(fe.vPrecise);
#endif
		r.low = 0;
		r.high = std::max(-fe.low, fe.high);
		//r.Check();
		return r;
	}
}

FORCE_INLINE bool Quadratic(EReal A, EReal B, EReal C, EReal *t0, EReal *t1);
FORCE_INLINE bool Quadratic(EReal A, EReal B, EReal C, EReal *t0, EReal *t1) {
	// Find quadratic discriminant
	double discrim = (double)B.v * (double)B.v - 4. * (double)A.v * (double)C.v;
	if (discrim < 0.) return false;
	double rootDiscrim = std::sqrt(discrim);

	EReal floatRootDiscrim(rootDiscrim, MachineEps * rootDiscrim);

	// Compute quadratic _t_ values
	EReal q;
	if ((real)B < 0)
		q = -.5 * (B - floatRootDiscrim);
	else
		q = -.5 * (B + floatRootDiscrim);
	*t0 = q / A;
	*t1 = C / q;
	if ((real)*t0 > (real)*t1) std::swap(*t0, *t1);
	return true;
}

GY_NAMESPACE_END