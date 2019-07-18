#pragma once

#include <string.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <functional>
#include "utils.h"
#include "sse.h"
#include "avx.h"

NAMESPACE_BEGIN

enum class IntrinsicSet {
	None,
	SSE,
	AVX
};

#if defined ISE_NONE
constexpr IntrinsicSet defaultInstructionSet = IntrinsicSet::None;
#elif defined ISE_SSE
constexpr IntrinsicSet defaultInstructionSet = IntrinsicSet::SSE;
#define USING_SIMD
#elif defined ISE_AVX
#define USING_SIMD
#endif

template<int dim, typename T, IntrinsicSet ISE, class Enable = void>
struct VectorBase {
	static constexpr int elements = dim;
	static constexpr bool simd = false;
	T d[dim];
};

template<typename T, IntrinsicSet ISE>
struct VectorBase<1, T, ISE, void> {
	static constexpr int elements = 1;
	static constexpr bool simd = false;
	union {
		T d[1];
		struct {
			T x;
		};
	};
};

template<typename T, IntrinsicSet ISE>
struct VectorBase<2, T, ISE, void> {
	static constexpr int elements = 2;
	static constexpr int simd = false;
	union {
		T d[2];
		struct {
			T x, y;
		};
	};
};


template<typename T, IntrinsicSet ISE>
struct VectorBase<3, T, ISE, typename std::enable_if_t<ISE < IntrinsicSet::SSE>> {
	static constexpr int elements = 3;
	static constexpr bool simd = false;
	union {
		T d[3];
		struct {
			T x, y, z;
		};
	};
};

template<IntrinsicSet ISE>
struct ALIGNED(16) VectorBase<3, float32, ISE, typename std::enable_if_t<ISE == IntrinsicSet::SSE> > {
	static constexpr int elements = 3;
	static constexpr bool simd = true;
	union {
		__m128 v;
		float32 d[4];
		struct {
			float32 x, y, z, _w;
		};
	};

	FORCE_INLINE VectorBase(float32 a = 0.f): v(_mm_set_ps1(a)) {}

	explicit VectorBase(__m128 v) : v(v) {}
};

template<typename T, IntrinsicSet ISE>
struct VectorBase<4, T, ISE, std::enable_if_t<ISE < IntrinsicSet::SSE>> {
	static constexpr int elements = 4;
	static constexpr bool simd = false;
	union {
		T d[4];
		struct {
			T x, y, z, w;
		};
	};
};

template<IntrinsicSet ISE>
struct ALIGNED(16) VectorBase < 4, float32, ISE, std::enable_if_t<ISE >= IntrinsicSet::SSE> > {
	static constexpr int elements = 4;
	static constexpr bool simd = true;
	union {
		__m128 v;
		float32 d[4];
		struct {
			float32 x, y, z, w;
		};
	};

	FORCE_INLINE VectorBase(float32 a = 0.f) : v(_mm_set_ps1(a)) {}

	explicit VectorBase(__m128 v) : v(v) {}
};


template<int dim_, typename T, IntrinsicSet ISE = defaultInstructionSet>
struct Vector : public VectorBase<dim_, T, ISE> {
	static constexpr int dim = dim_;
	static constexpr int elements = dim_;

	template <int dim__, typename T_, IntrinsicSet ISE_>
	static constexpr bool SIMD_4_32F = (dim__ == 3 || dim__ == 4) &&
		std::is_same<T_, float32>::value && ISE_ >= IntrinsicSet::SSE;

	template <int dim__, typename T_, IntrinsicSet ISE_>
	static constexpr bool SIMD_4_64F = (dim__ == 3 || dim__ == 4) &&
		std::is_same<T_, float32>::value && ISE_ >= IntrinsicSet::SSE;

	template <int dim__, typename T_, IntrinsicSet ISE_>
	static constexpr bool NOT_SIMD_4_32F = !SIMD_4_32F<dim__, T_, ISE_>;

	template <int dim__, typename T_, IntrinsicSet ISE_>
	static constexpr bool SIMD_NONE = !SIMD_4_32F<dim__, T_, ISE_>;

	template<int dim__, typename T_, IntrinsicSet ISE_>
	static constexpr bool SSE_4_32F = (dim__ == 3 || dim__ == 4) &&
		std::is_same<T_, float32>::value && ISE_ == IntrinsicSet::SSE;
	
	using type = T;

	FORCE_INLINE Vector() {
		for (int i = 0; i < dim; ++i) {
			this->d[i] = T(0);
		}
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SSE_4_32F<dim__, T_, ISE_>, int> = 0>
	explicit FORCE_INLINE Vector() : VectorBase(_mm_set_ps1(0.f)) { }

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SSE_4_32F<dim__, T_, ISE_>, int> = 0 >
	explicit FORCE_INLINE Vector(__m128 v) : VectorBase<dim__, T_, ISE_>(v) { }

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		std::enable_if_t<SSE_4_32F<dim__, T_, ISE_>, int> = 0>
	explicit FORCE_INLINE Vector(real x, real y, real z, real w = 0.f) :
		VectorBase<dim__, T_, ISE_>(_mm_set_ps(w, z, y, x)) { }


	FORCE_INLINE Vector(T *a) {
		memcpy(&(this->d[0]), a, sizeof(T) * dim);
	}

	
	template<typename T_>
	FORCE_INLINE Vector(const std::vector<T_> &a) {
		for (int i = 0; i < dim; ++i) {
			this->d[i] = a[i];
		}
	}

	FORCE_INLINE Vector(const Vector &a) {
		for (int i = 0; i < dim; ++i) {
			this->d[i] = a.d[i];
		}
	}

	FORCE_INLINE Vector(Vector &&a) {
		for (int i = 0; i < dim; ++i) {
			this->d[i] = std::move(a.d[i]);
		}
	}

	FORCE_INLINE Vector(T a) {
		for (int i = 0; i < dim; ++i) {
			this->d[i] = a;
		}
	}

	template<int dim__ = dim, typename T_ = T, typename std::enable_if_t<dim__ == 1, int> = 0>
	FORCE_INLINE Vector(T x) {
		this->d[0] = x;
	}

	template<int dim__ = dim, typename T_ = T, typename std::enable_if_t<dim__ == 2, int> = 0>
	FORCE_INLINE Vector(T x, T y) {
		this->d[0] = x;
		this->d[1] = y;
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE, 
		typename std::enable_if_t<dim__ == 3 && ((ISE_ == IntrinsicSet::None) || (ISE_ >= IntrinsicSet::SSE && !std::is_same<T_, float32>::value))
			, int> = 0>
	FORCE_INLINE Vector(T x, T y, T z) {
		this->d[0] = x;
		this->d[1] = y;
		this->d[2] = z;
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<dim__ == 4 && ((ISE_ == IntrinsicSet::None) || (ISE_ >= IntrinsicSet::SSE && !std::is_same<T_, float32>::value))
		, int> = 0>
	FORCE_INLINE Vector(T x, T y, T z, T w) {
		this->d[0] = x;
		this->d[1] = y;
		this->d[2] = z;
		this->d[3] = w;
	}

	FORCE_INLINE T& operator[](int i) {
		return this->d[i];
	}

	FORCE_INLINE const T& operator[](int i) const {
		return this->d[i];
	}

	FORCE_INLINE T& operator()(int i) {
		return this->d[i];
	}

	FORCE_INLINE const T& operator()(int i) const {
		return this->d[i];
	}

	FORCE_INLINE Vector& operator=(const Vector &a) {
		//memcpy(&(this->d[0]), &(a.d[0]), sizeof(T) * dim);
		memcpy(this, &a, sizeof(*this));
		return *this;
	}

	FORCE_INLINE Vector& operator=(Vector &&a) {
		//memcpy(&(this->d[0]), &(a.d[0]), sizeof(T) * dim);
		memcpy(this, &a, sizeof(*this));
		return *this;
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SIMD_NONE<dim__, T_, ISE_>, int> = 0 >
	FORCE_INLINE Vector operator+(const Vector &a) const {
		Vector ret;
		for (int i = 0; i < dim; ++i) {
			ret.d[i] = this->d[i] + a.d[i];
		}
		return ret;
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SSE_4_32F<dim__, T_, ISE_>, int> = 0>
	FORCE_INLINE Vector operator+(const Vector &a) const {
		return Vector(_mm_add_ps(this->v, a.v));
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SIMD_NONE<dim__, T_, ISE_>, int> = 0 >
	FORCE_INLINE Vector operator-(const Vector &a) const {
		Vector ret;
		for (int i = 0; i < dim; ++i) {
			ret.d[i] = this->d[i] - a.d[i];
		}
		return ret;
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SSE_4_32F<dim__, T_, ISE_>, int> = 0>
	FORCE_INLINE Vector operator-(const Vector &a) const {
		return Vector(_mm_sub_ps(this->v, a.v));
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SIMD_NONE<dim__, T_, ISE_>, int> = 0 >
	FORCE_INLINE Vector operator-() const {
		Vector ret;
		for (int i = 0; i < dim; ++i) {
			ret.d[i] = -this->d[i];
		}
		return ret;
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SSE_4_32F<dim__, T_, ISE_>, int> = 0 >
	FORCE_INLINE Vector operator-() const {
		return Vector(_mm_sub_ps(_mm_set_ps1(0.f), this->v));
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SIMD_NONE<dim__, T_, ISE_>, int> = 0 >
	FORCE_INLINE Vector operator*(const Vector &a) const {
		Vector ret;
		for (int i = 0; i < dim; ++i) {
			ret.d[i] = this->d[i] * a.d[i];
		}
		return ret;
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SSE_4_32F<dim__, T_, ISE_>, int> = 0>
	FORCE_INLINE Vector operator*(const Vector &a) const {
		return Vector(_mm_mul_ps(this->v, a.v));
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SIMD_NONE<dim__, T_, ISE_>, int> = 0 >
	FORCE_INLINE Vector operator*(real a) const {
		Vector ret;
		for (int i = 0; i < dim; ++i) {
			ret.d[i] = this->d[i] * a;
		}
		return ret;
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SSE_4_32F<dim__, T_, ISE_>, int> = 0 >
	FORCE_INLINE Vector operator*(float32 a) const {
		return Vector(_mm_mul_ps(this->v, _mm_set_ps1(a)));
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SIMD_NONE<dim__, T_, ISE_>, int> = 0 >
	FORCE_INLINE Vector operator/(const Vector &a) const {
		Vector ret;
		for (int i = 0; i < dim; ++i) {
			ret.d[i] = this->d[i] / a.d[i];
		}
		return ret;
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SSE_4_32F<dim__, T_, ISE_>, int> = 0 >
	FORCE_INLINE Vector operator/(const Vector &a) const {
		return Vector(_mm_div_ps(this->v, a.v));
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SIMD_NONE<dim__, T_, ISE_>, int> = 0 >
	FORCE_INLINE Vector operator/(real a) const {
		Vector ret;
		if (a == 0) {
			return ret;
		}
		for (int i = 0; i < dim; ++i) {
			ret.d[i] = this->d[i] / a;
		}
		return ret;
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SSE_4_32F<dim__, T_, ISE_>, int> = 0 >
	FORCE_INLINE Vector operator/(float32 a) const {
		if(a == 0) return Vector(_mm_set_ps1(0.f));
		return Vector(_mm_div_ps(this->v, _mm_set_ps1(a)));
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SIMD_NONE<dim__, T_, ISE_>, int> = 0 >
	FORCE_INLINE Vector& operator+=(const Vector &a) {
		for (int i = 0; i < dim; ++i) {
			this->d[i] += a.d[i];
		}
		return *this;
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SSE_4_32F<dim__, T_, ISE_>, int> = 0>
	FORCE_INLINE Vector& operator+=(const Vector &a) {
		this->v = _mm_add_ps(this->v, a.v);
		return *this;
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SIMD_NONE<dim__, T_, ISE_>, int> = 0 >
	FORCE_INLINE Vector& operator-=(const Vector &a) {
		for (int i = 0; i < dim; ++i) {
			this->d[i] -= a.d[i];
		}
		return *this;
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SSE_4_32F<dim__, T_, ISE_>, int> = 0>
		FORCE_INLINE Vector& operator-=(const Vector &a) {
		this->v = _mm_sub_ps(this->v, a.v);
		return *this;
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SIMD_NONE<dim__, T_, ISE_>, int> = 0 >
	FORCE_INLINE Vector& operator*=(const Vector &a) {
		for (int i = 0; i < dim; ++i) {
			this->d[i] *= a.d[i];
		}
		return *this;
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SSE_4_32F<dim__, T_, ISE_>, int> = 0>
	FORCE_INLINE Vector& operator*=(const Vector &a) {
		this->v = _mm_mul_ps(this->v, a.v);
		return *this;
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SIMD_NONE<dim__, T_, ISE_>, int> = 0 >
	FORCE_INLINE Vector& operator*=(real a) {
		for (int i = 0; i < dim; ++i) {
			this->d[i] *= a;
		}
		return *this;
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SSE_4_32F<dim__, T_, ISE_>, int> = 0>
	FORCE_INLINE Vector& operator*=(real a) {
		this->v = _mm_mul_ps(this->v, _mm_set_ps1(a));
		return *this;
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SIMD_NONE<dim__, T_, ISE_>, int> = 0 >
	FORCE_INLINE Vector& operator/=(const Vector &a) {
		for (int i = 0; i < dim; ++i) {
			this->d[i] /= a.d[i];
		}
		return *this;
	}


	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SSE_4_32F<dim__, T_, ISE_>, int> = 0>
	FORCE_INLINE Vector& operator/=(const Vector &a) {
		this->v = _mm_div_ps(this->v, a.v);
		return *this;
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SIMD_NONE<dim__, T_, ISE_>, int> = 0 >
	FORCE_INLINE Vector& operator/=(real a) {
		for (int i = 0; i < dim; ++i) {
			this->d[i] /= a;
		}
		return *this;
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SSE_4_32F<dim__, T_, ISE_>, int> = 0 >
	FORCE_INLINE Vector& operator/=(real a) {
		this->v = _mm_div_ps(this->v, _mm_set_ps1(a));
		return *this;
	}

	FORCE_INLINE bool operator==(const Vector &a) const {
		for (int i = 0; i < dim; ++i) {
			if (this->d[i] != a.d[i]) {
				return false;
			}
		}
		return true;
	}

	FORCE_INLINE bool operator!=(const Vector &a) const {
		for (int i = 0; i < dim; ++i) {
			if (this->d[i] != a.d[i]) {
				return true;
			}
		}
		return false;
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SIMD_NONE<dim__, T_, ISE_>, int> = 0 >
	FORCE_INLINE T Dot(const Vector &a) const {
		T ret = T(0);
		for (int i = 0; i < dim; ++i) {
			ret += this->d[i] * a.d[i];
		}
		return ret;
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SSE_4_32F<dim__, T_, ISE_> && dim__ == 3, int> = 0>
	FORCE_INLINE real Dot(const Vector &a) const {
		return _mm_cvtss_f32(_mm_dp_ps(this->v, a.v, 0x71));
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SSE_4_32F<dim__, T_, ISE_> && dim__ == 4, int> = 0>
	FORCE_INLINE real Dot(const Vector &a) const {
		return _mm_cvtss_f32(_mm_dp_ps(this->v, a.v, 0xf1));
	}

	FORCE_INLINE real Length2() const {
		return this->Dot(*this);
	}

	FORCE_INLINE real Length() const {
		return std::sqrt(Length2());
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SIMD_NONE<dim__, T_, ISE_>, int> = 0 >
	FORCE_INLINE Vector Norm() const {
		real invLen = (real)(1.0) / Length();
		Vector ret(*this);
		for (int i = 0; i < dim; ++i) {
			ret.d[i] *= invLen;
		}
		return ret;
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SSE_4_32F<dim__, T_, ISE_>, int> = 0>
	FORCE_INLINE Vector Norm() const {
		return Vector(_mm_div_ps(this->v, _mm_set_ps1(this->Length())));
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SIMD_NONE<dim__, T_, ISE_>, int> = 0 >
	FORCE_INLINE void Normalize() {
		real invLen = (real)(1.0) / Length();
		for (int i = 0; i < dim; ++i) {
			this->d[i] *= invLen;
		}
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SSE_4_32F<dim__, T_, ISE_>, int> = 0 >
	FORCE_INLINE void Normalize() {
		this->v = _mm_div_ps(this->v, _mm_set_ps1(this->Length()));
	}

	FORCE_INLINE T MaxVal() const {
		T ret = this->d[0];
		for (int i = 1; i < dim; ++i) {
			ret = std::max(ret, this->d[i]);
		}
		return ret;
	}

	FORCE_INLINE T MinVal() const {
		T ret = this->d[0];
		for (int i = 1; i < dim; ++i) {
			ret = std::min(ret, this->d[i]);
		}
		return ret;
	}

	template<int dim__ = dim, typename T, typename std::enable_if_t<dim__ == 3, int> = 0>
	FORCE_INLINE Vector Cross(const Vector<dim__, T> &a) const {
		return Vector<dim__, T>(this->y * a.z - this->z * a.y, this->z * a.x - this->x * a.z, this->x * a.y - this->y * a.x);
	}

	FORCE_INLINE real Y()  const {
		const real YWeight[3] = { (real)0.212671, (real)0.715160, (real)0.072169f };
		return YWeight[0] * this->x + YWeight[1] * this->y + YWeight[2] * this->z;
	}

	template <int a,
		int b,
		int c,
		int d,
		int dim__ = dim,
		typename T_ = T,
		IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<!SIMD_4_32F<dim__, T_, ISE_>, int> = 0>
	FORCE_INLINE Vector Permute() const {
		return Vector(this->d[a], this->d[b], this->d[c], this->d[d]);
	}

	template <int a,
		int b,
		int c,
		int d,
		int dim__ = dim,
		typename T_ = T,
		IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SIMD_4_32F<dim__, T_, ISE_>, int> = 0>
	FORCE_INLINE Vector Permute() const {
		return Vector(_mm_permute_ps(this->v, _MM_SHUFFLE(a, b, c, d)));
	}

	template <int a, int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE>
	FORCE_INLINE Vector Broadcast() const {
		return Permute<a, a, a, a>();
	}
};

template<int dim, typename T, IntrinsicSet ISE>
FORCE_INLINE std::ostream& operator<<(std::ostream &os, const Vector<dim, T, ISE> &a) {
	os << "[ ";
	for (int i = 0; i < dim; ++i) {
		if (i < dim - 1) {
			os << a.d[i] << ", ";
		}
		else {
			os << a.d[i] << " ";
		}

	}
	os << "]";
	return os;
}

template<int dim, typename T, IntrinsicSet ISE>
FORCE_INLINE Vector<dim, T, ISE> operator*(real a, const Vector<dim, T, ISE> &b) {
	return b * a;
}


template<int dim, typename T, IntrinsicSet ISE>
FORCE_INLINE real Dot(const Vector<dim, T, ISE> &a, const Vector<dim, T, ISE> &b) {
	return a.Dot(b);
}

template<typename T, IntrinsicSet ISE>
FORCE_INLINE Vector<3, T, ISE> Cross(const Vector<3, T, ISE> &a, const Vector<3, T, ISE> &b) {
	return Vector<3, T, ISE>(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

template<int dim, typename T, IntrinsicSet ISE>
FORCE_INLINE void Normalize(Vector<dim, T, ISE> &a) {
	a.Normalize();
}

template<int dim, typename T, IntrinsicSet ISE>
FORCE_INLINE real Distance2(const Vector<dim, T, ISE> &a, const Vector<dim, T, ISE> &b) {
	return (b - a).Length2();
}

template<int dim, typename T, IntrinsicSet ISE>
FORCE_INLINE real Distance(const Vector<dim, T, ISE> &a, const Vector<dim, T, ISE> &b) {
	return (b - a).Length();
}

template<int dim, typename T, IntrinsicSet ISE>
FORCE_INLINE Vector<dim, T, ISE> Abs(const Vector<dim, T, ISE> &a) {
	Vector<dim, T, ISE> ret;
	for(int i = 0; i < dim; ++i) {
		ret.d[i] = std::abs(a.d[i]);
	}
	return ret;
}

using Vector2 = Vector<2, real, defaultInstructionSet>;
using Vector2i = Vector<2, int, defaultInstructionSet>;
using Vector2f = Vector<2, float, defaultInstructionSet>;
using Vector2d = Vector<2, double, defaultInstructionSet>;

using Vector3 = Vector<3, real, defaultInstructionSet>;
using Vector3i = Vector<3, int, defaultInstructionSet>;
using Vector3f = Vector<3, float, defaultInstructionSet>;
using Vector3d = Vector<3, double, defaultInstructionSet>;

using Vector4 = Vector<4, real, defaultInstructionSet>;
using Vector4i = Vector<4, int, defaultInstructionSet>;
using Vector4f = Vector<4, float, defaultInstructionSet>;
using Vector4d = Vector<4, double, defaultInstructionSet>;

using Vec3 = Vector3;
using Vec2 = Vector2;


// FMA: a * b + c
template <typename T>
FORCE_INLINE typename std::enable_if<T::simd && (defaultInstructionSet >= IntrinsicSet::AVX), T>::type
fused_mul_add(const T &a, const T &b, const T &c) {
	return T(_mm_fmadd_ps(a, b, c));
}

template <typename T>
FORCE_INLINE typename std::enable_if<!T::simd || (defaultInstructionSet < IntrinsicSet::AVX), T>::type
fused_mul_add(const T &a, const T &b, const T &c) {
	return a * b + c;
}


template<int dim_, typename T, IntrinsicSet ISE>
struct Matrix {
	static constexpr int dim = dim_;
	using type = T;
	Vector<dim, T, ISE> d[dim];

	template<int dim__, typename T_, IntrinsicSet ISE_>
	static constexpr bool SIMD_4_32F =
		(dim__ == 3 || dim__ == 4) && std::is_same<T_, float32>::value &&
		(ISE_ >= IntrinsicSet::SSE);


	template <
		typename F,
		std::enable_if_t<std::is_convertible<
		F,
		std::function<Vector<dim_, T, ISE>(int)>>::value,
		int> = 0>
	FORCE_INLINE Matrix& Set(const F &f) {
		for (int i = 0; i < dim; i++)
			this->d[i] = f(i);
		return *this;
	}

	FORCE_INLINE Matrix() {
		for (int i = 0; i < dim; ++i) {
			d[i] = Vector<dim, T, ISE>();
			d[i][i] = 1;
		}
	}

	explicit FORCE_INLINE Matrix(const Matrix &m) {
		*this = m;
	}

	FORCE_INLINE Matrix(T a) {
		for (int i = 0; i < dim; ++i) {
			d[i][i] = a;
		}
	}

	FORCE_INLINE Matrix(const Vector<dim, T, ISE> &v) {
		for (int i = 0; i < dim; i++)
			this->d[i][i] = v[i];
	}

	template<int dim__ = dim, typename T_, 
		typename std::enable_if_t<dim__ == 2> = 0>
	FORCE_INLINE Matrix(const Vector<dim__, T_> &v0, const Vector<dim__, T_> &v1) {
		this->d[0] = v0;
		this->d[1] = v1;
	}

	template<int dim__ = dim, typename T_,
		typename std::enable_if_t<dim__ == 3> = 0>
	FORCE_INLINE Matrix(const Vector<dim__, T_> &v0, const Vector<dim__, T_> &v1, const Vector<dim__, T_> &v2) {
		this->d[0] = v0;
		this->d[1] = v1;
		this->d[2] = v2;
	}

	template<int dim__ = dim, typename T_,
		typename std::enable_if_t<dim__ == 4> = 0>
	FORCE_INLINE Matrix(const Vector<dim__, T_> &v0, const Vector<dim__, T_> &v1, 
		const Vector<dim__, T_> &v2, const Vector<dim__, T_> &v3) {
		this->d[0] = v0;
		this->d[1] = v1;
		this->d[2] = v2;
		this->d[3] = v3;
	}

	template <typename F, 
		typename std::enable_if_t<
		std::is_convertible<F, std::function<Vector<dim, T, ISE>(int)>>::value,
		int> = 0>
	explicit FORCE_INLINE Matrix(const F &f) {
		for (int i = 0; i < dim; i++)
			this->d[i] = f(i);
	}

	template<int dim__ = dim, typename std::enable_if_t<dim__ == 2, int> = 0>
	FORCE_INLINE Matrix(const Vector<dim, T> &v0, const Vector<dim, T> &v1) {
		d[0] = v0;
		d[1] = v1;
	}

	template<int dim__ = dim, typename std::enable_if_t<dim__ == 2, int> = 0>
	FORCE_INLINE Matrix(T m00, T m01,
						T m10, T m11) {
		d[0][0] = m00; d[1][0] = m01;
		d[0][1] = m10; d[1][1] = m11;
	}

	template<int dim__ = dim, typename std::enable_if_t<dim_ == 3, int> = 0>
	FORCE_INLINE Matrix(const Vector<dim, T> &v0, const Vector<dim, T> &v1, const Vector<dim, T> &v2) {
		d[0] = v0;
		d[1] = v1;
		d[2] = v2;
	}

	template<int dim__ = dim, typename = std::enable_if_t<dim__ == 3, int> = 0>
	FORCE_INLINE Matrix(T m00, T m01, T m02,
						T m10, T m11, T m12,
						T m20, T m21, T m22) {
		d[0][0] = m00; d[1][0] = m01; d[2][0] = m02;
		d[0][1] = m10; d[1][1] = m11; d[2][1] = m12;
		d[0][2] = m20; d[1][2] = m21; d[2][2] = m22;
	}

	template<int dim__ = dim, typename std::enable_if_t<dim__ == 4, int> = 0>
	FORCE_INLINE Matrix(const Vector<dim, T> &v0, const Vector<dim, T> &v1, 
		const Vector<dim, T> &v2, const Vector<dim, T> &v3) {
		d[0] = v0;
		d[1] = v1;
		d[2] = v2;
		d[3] = v3;
	}

	template<int dim__ = dim, typename std::enable_if_t<dim__ == 3, int> = 0>
	FORCE_INLINE Matrix(T m00, T m01, T m02, T m03,
						T m10, T m11, T m12, T m13,
						T m20, T m21, T m22, T m23,
						T m30, T m31, T m32, T m33) {
		d[0][0] = m00; d[1][0] = m01; d[2][0] = m02; d[3][0] = m03;
		d[0][1] = m10; d[1][1] = m11; d[2][1] = m12; d[3][1] = m13;
		d[2][2] = m20; d[1][2] = m21; d[2][2] = m22; d[3][2] = m23;
		d[0][3] = m30; d[1][3] = m31; d[2][3] = m32; d[3][3] = m33;
	}

	FORCE_INLINE Matrix& operator=(const Matrix &m) {
		for (int i = 0; i < dim; ++i) {
			this->d[i] = m.d[i];
		}
		return *this;
	}

	FORCE_INLINE Matrix& operator=(Matrix &&m) {
		for (int i = 0; i < dim; ++i) {
			this->d[i] = m.d[i];
		}
		return *this;
	}

	FORCE_INLINE Vector<dim, T>& operator[](int i) {
		return d[i];
	}

	FORCE_INLINE const T& operator()(int i, int j) const {
		return d[j][i];
	}

	FORCE_INLINE T& operator()(int i, int j) {
		return d[j][i];
	}

	FORCE_INLINE Matrix operator+(const Matrix &a) const {
		return Matrix([=](int i) {this->d[i] + a.d[i]; });
	}

	FORCE_INLINE Matrix& operator+=(const Matrix &a) {
		return this->Set([&](int i) {return this->d[i] + a[i]; });
	}

	FORCE_INLINE Matrix operator-(const Matrix &a) const {
		return Matrix([=](int i) {this->d[i] - a.d[i]; });
	}

	FORCE_INLINE Matrix& operator-=(const Matrix &a) {
		return this->Set([&](int i) {return this->d[i] - a[i]; });
	}

	FORCE_INLINE bool operator==(const Matrix &a) const {
		for (int i = 0; i < dim; i++)
			for (int j = 0; j < dim; j++)
				if (d[i][j] != a[i][j])
					return false;
		return true;
	}

	FORCE_INLINE bool operator!=(const Matrix &o) const {
		for (int i = 0; i < dim; i++)
			for (int j = 0; j < dim; j++)
				if (d[i][j] != o[i][j])
					return true;
		return false;
	}

	template<int dim__ = dim, typename T_ = T, IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<!SIMD_4_32F<dim__, T_, ISE_>, int> = 0>
	FORCE_INLINE Vector<dim, T, ISE> operator*(const Vector<dim, T, ISE> &a) const {
		Vector<dim, T, ISE> ret(d[0] * a.d[0]);
		for (int i = 1; i < dim; ++i) {
			ret += d[i] * a.d[i];
		}
		return ret;
	}

	// Matrix3
	template <int dim__ = dim,
		typename T_ = T,
		IntrinsicSet ISE_ = ISE,
		typename std::enable_if_t<SIMD_4_32F<dim__, T_, ISE_> && dim__ == 3,
		int> = 0>
	FORCE_INLINE Vector<dim, T, ISE> operator*(const Vector<dim, T, ISE> &a) const {
		Vector<dim, T, ISE> ret = a.template broadcast<2>() * d[2];
		ret = fused_mul_add(d[1], a.template broadcast<1>(), ret);
		ret = fused_mul_add(d[0], a.template broadcast<0>(), ret);
		return ret;
	}

	FORCE_INLINE Matrix operator*(const Matrix &a) const {
		return Matrix([=](int i) {(*this) * a[i]; });
	}

	FORCE_INLINE Matrix operator*(real a) const {
		return Matrix([=](int i) {(*this)[i] * a; });
	}

	FORCE_INLINE Matrix& operator*=(const Matrix &a) {
		Matrix mat([&](int i) {(*this) * a[i]; });
		*this = std::move(mat);
		return *this;
	}

	FORCE_INLINE Matrix& operator*=(real a) {
		for (int i = 0; i < dim; ++i) {
			d[i] *= a;
		}
		return *this;
	}

	FORCE_INLINE Matrix operator/(real a) const {
		return Matrix([=](int i) {(*this)[i] / a; });
	}

	FORCE_INLINE Matrix& operator/=(real a) {
		for (int i = 0; i < dim; ++i) {
			d[i] /= a;
		}
		return *this;
	}

	FORCE_INLINE T FrobeniusNorm2() const {
		T sum = d[0].Length2();
		for (int i = 1; i < dim; i++) {
			sum += d[i].Length2();
		}
		return sum;
	}

	FORCE_INLINE auto FrobeniusNorm() const {
		return std::sqrt(FrobeniusNorm2());
	}

	FORCE_INLINE Matrix Transpose() const {
		Matrix ret;
		for (int i = 0; i < dim; i++) {
			for (int j = 0; j < dim; j++) {
				ret[i][j] = d[j][i];
			}
		}
		return ret;
	}

	FORCE_INLINE static Matrix Identidy() {
		return Matrix(1.f);
	}

	FORCE_INLINE static Matrix Zeros() {
		return Matrix();
	}
};


using Matrix2 = Matrix<2, real, defaultInstructionSet>;
using Matrix3 = Matrix<3, real, defaultInstructionSet>;
using Matrix4 = Matrix<4, real, defaultInstructionSet>;

using Matrix2f = Matrix<2, float32, defaultInstructionSet>;
using Matrix3f = Matrix<3, float32, defaultInstructionSet>;
using Matrix4f = Matrix<4, float32, defaultInstructionSet>;

using Matrix2d = Matrix<2, float64, defaultInstructionSet>;
using Matrix3d = Matrix<3, float64, defaultInstructionSet>;
using Matrix4d = Matrix<4, float64, defaultInstructionSet>;

NAMESPACE_END