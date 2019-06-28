#pragma once

#include "utils.h"
#include <string.h>
#include <iostream>
#include <algorithm>
#include <vector>

NAMESPACE_BEGIN

template<int dim, typename T>
struct VectorBase {
	static constexpr int elements = dim;
	T d[dim];
};

template<typename T>
struct VectorBase<1, T> {
	static constexpr int elements = 1;
	union {
		T d[1];
		struct {
			T x;
		};
	};
};

template<typename T>
struct VectorBase<2, T> {
	static constexpr int elements = 2;
	union {
		T d[2];
		struct {
			T x, y;
		};
	};
};

template<typename T>
struct VectorBase<3, T> {
	static constexpr int elements = 3;
	union {
		T d[3];
		struct {
			T x, y, z;
		};
	};
};

template<typename T>
struct VectorBase<4, T> {
	static constexpr int elements = 4;
	union {
		T d[4];
		struct {
			T x, y, z, w;
		};
	};
};


template<int dim_, typename T>
struct Vector : public VectorBase<dim_, T> {
	static constexpr int dim = dim_;
	static constexpr int elements = dim_;

	using type = T;

	FORCE_INLINE Vector() {
		for (int i = 0; i < dim; ++i) {
			this->d[i] = T(0);
		}
	}

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

	template<int dim__ = dim, typename T_ = T, typename = std::enable_if_t<dim == 1, int>>
	FORCE_INLINE Vector(T x) {
		this->d[0] = x;
	}

	template<int dim__ = dim, typename T_ = T, typename = std::enable_if_t<dim == 2, int>>
	FORCE_INLINE Vector(T x, T y) {
		this->d[0] = x;
		this->d[1] = y;
	}

	template<int dim__ = dim, typename T_ = T, typename = std::enable_if_t<dim == 3, int>>
	FORCE_INLINE Vector(T x, T y, T z) {
		this->d[0] = x;
		this->d[1] = y;
		this->d[2] = z;
	}

	template<int dim__ = dim, typename T_ = T, typename = std::enable_if_t<dim == 4, int>>
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
		memcpy(&(this->d[0]), &(a.d[0]), sizeof(T) * dim);
		return *this;
	}

	FORCE_INLINE Vector& operator=(Vector &&a) {
		memcpy(&(this->d[0]), &(a.d[0]), sizeof(T) * dim);
		return *this;
	}

	FORCE_INLINE Vector operator+(const Vector &a) const {
		Vector ret;
		for (int i = 0; i < dim; ++i) {
			ret.d[i] = this->d[i] + a.d[i];
		}
		return ret;
	}

	FORCE_INLINE Vector operator-(const Vector &a) const {
		Vector ret;
		for (int i = 0; i < dim; ++i) {
			ret.d[i] = this->d[i] - a.d[i];
		}
		return ret;
	}

	FORCE_INLINE Vector operator-() const {
		Vector ret;
		for (int i = 0; i < dim; ++i) {
			ret.d[i] = -this->d[i];
		}
		return ret;
	}

	FORCE_INLINE Vector operator*(const Vector &a) const {
		Vector ret;
		for (int i = 0; i < dim; ++i) {
			ret.d[i] = this->d[i] * a.d[i];
		}
		return ret;
	}

	FORCE_INLINE Vector operator*(real a) const {
		Vector ret;
		for (int i = 0; i < dim; ++i) {
			ret.d[i] = this->d[i] * a;
		}
		return ret;
	}

	FORCE_INLINE Vector operator/(const Vector &a) const {
		Vector ret;
		for (int i = 0; i < dim; ++i) {
			ret.d[i] = this->d[i] / a.d[i];
		}
		return ret;
	}

	FORCE_INLINE Vector operator/(real a) const {
		Vector ret;
		for (int i = 0; i < dim; ++i) {
			ret.d[i] = this->d[i] / a;
		}
		return ret;
	}

	FORCE_INLINE Vector& operator+=(const Vector &a) {
		for (int i = 0; i < dim; ++i) {
			this->d[i] += a.d[i];
		}
		return *this;
	}

	FORCE_INLINE Vector& operator-=(const Vector &a) {
		for (int i = 0; i < dim; ++i) {
			this->d[i] -= a.d[i];
		}
		return *this;
	}

	FORCE_INLINE Vector& operator*=(const Vector &a) {
		for (int i = 0; i < dim; ++i) {
			this->d[i] *= a.d[i];
		}
		return *this;
	}

	FORCE_INLINE Vector& operator*=(real a) {
		for (int i = 0; i < dim; ++i) {
			this->d[i] *= a;
		}
		return *this;
	}

	FORCE_INLINE Vector& operator/=(const Vector &a) {
		for (int i = 0; i < dim; ++i) {
			this->d[i] /= a.d[i];
		}
		return *this;
	}

	FORCE_INLINE Vector& operator/=(real a) {
		for (int i = 0; i < dim; ++i) {
			this->d[i] /= a;
		}
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

	FORCE_INLINE T Dot(const Vector &a) const {
		T ret = T(0);
		for (int i = 0; i < dim; ++i) {
			ret += this->d[i] * a.d[i];
		}
		return ret;
	}

	FORCE_INLINE real Length2() const {
		real ret = real(0);
		for (int i = 0; i < dim; ++i) {
			ret += this->d[i] * this->d[i];
		}
		return ret;
	}

	FORCE_INLINE real Length() const {
		return std::sqrt(Length2());
	}

	FORCE_INLINE Vector Normal() const {
		real invLen = (real)(1.0) / Length();
		Vector ret(*this);
		for (int i = 0; i < dim; ++i) {
			ret.d[i] *= invLen;
		}
		return ret;
	}

	FORCE_INLINE void Normalize() {
		real invLen = (real)(1.0) / Length();
		for (int i = 0; i < dim; ++i) {
			this->d[i] *= invLen;
		}
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

	template<int dim__ = dim, typename T, typename = std::enable_if_t<dim__ == 3, int>>
	FORCE_INLINE Vector Cross(const Vector<dim__, T> &a) const {
		return Vector<dim__, T>(this->y * a.z - this->z * a.y, this->z * a.x - this->x * a.z, this->x * a.y - this->y * a.x);
	}

	FORCE_INLINE real Y()  const {
		const real YWeight[3] = { (real)0.212671, (real)0.715160, (real)0.072169f };
		return YWeight[0] * this->x + YWeight[1] * this->y + YWeight[2] * this->z;
	}
};

template<int dim, typename T>
FORCE_INLINE std::ostream& operator<<(std::ostream &os, const Vector<dim, T> &a) {
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

template<int dim, typename T>
FORCE_INLINE Vector<dim, T> operator*(real a, const Vector<dim, T> &b) {
	Vector<dim, T> ret;
	for (int i = 0; i < dim; ++i) {
		ret.d[i] = a * b.d[i];
	}
	return ret;
}

template<int dim, typename T>
FORCE_INLINE T Dot(const Vector<dim, T> &a, const Vector<dim, T> &b) {
	T ret = T(0);
	for (int i = 0; i < dim; ++i) {
		ret += a.d[i] * b.d[i];
	}
	return ret;
}

template<typename T>
FORCE_INLINE Vector<3, T> Cross(const Vector<3, T> &a, const Vector<3, T> &b) {
	return Vector<3, T>(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
}

template<int dim, typename T>
FORCE_INLINE void Normalize(Vector<dim, T> &a) {
	a.Normalize();
}

template<int dim, typename T>
FORCE_INLINE real Distance(const Vector<dim, T> &a, const Vector<dim, T> &b) {
	return (b - a).Length();
}

using Vector2 = Vector<2, real>;
using Vector2i = Vector<2, int>;
using Vector2f = Vector<2, float>;
using Vector2d = Vector<2, double>;

using Vector3 = Vector<3, real>;
using Vector3i = Vector<3, int>;
using Vector3f = Vector<3, float>;
using Vector3d = Vector<3, double>;

using Vector4 = Vector<4, real>;
using Vector4i = Vector<4, int>;
using Vector4f = Vector<4, float>;
using Vector4d = Vector<4, double>;



struct Vec {
	real x, y, z; // vector: position, also color (r,g,b)
	Vec(real x_ = 0, real y_ = 0, real z_ = 0) { x = x_; y = y_; z = z_; }
	inline Vec operator+(const Vec &b) const { return Vec(x + b.x, y + b.y, z + b.z); }
	inline Vec operator-(const Vec &b) const { return Vec(x - b.x, y - b.y, z - b.z); }
	inline Vec operator+(real b) const { return Vec(x + b, y + b, z + b); }
	inline Vec operator-(real b) const { return Vec(x - b, y - b, z - b); }
	inline Vec operator*(real b) const { return Vec(x * b, y * b, z * b); }
	inline Vec operator*(const Vec &b) const { return Vec(x * b.x, y * b.y, z * b.z); }
	inline Vec operator/(real b) const { if (b == 0) return Vec(); else return Vec(x / b, y / b, z / b); }
	inline bool operator==(const Vec &b) const { return x == b.x && y == b.y && z == b.z; }
	inline bool operator!=(const Vec &b) const { return x != b.x || y != b.y || z != b.z; }
	inline Vec mul(const Vec &b) const { return Vec(x * b.x, y * b.y, z * b.z); }
	inline Vec norm() { return (*this) * (1.0 / sqrt(x * x + y * y + z * z)); }
	inline void normalize() {
		real invLength = 1.0 / std::sqrt(x * x + y * y + z * z);
		this->x *= invLength;
		this->y *= invLength;
		this->z *= invLength;
	}
	inline real dot(const Vec &b) const { return x * b.x + y * b.y + z * b.z; }
	inline real length() const { return std::sqrt(x * x + y * y + z * z); }
	inline real length2() const { return x * x + y * y + z * z; }
	Vec operator%(Vec&b) const { return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }
	real& operator[](int i) { return i == 0 ? x : i == 1 ? y : z; }
	real maxValue() const {
		return std::max(x, std::max(y, z));
	}
	real Y() const {
		const real YWeight[3] = { 0.212671f, 0.715160f, 0.072169f };
		return YWeight[0] * x + YWeight[1] * y + YWeight[2] * z;
	}
};

Vec operator*(real a, Vec b);

std::ostream& operator<<(std::ostream &os, const Vec &v);




template<int dim_, typename T>
struct Matrix {
	static constexpr int dim = dim_;
	using type = T;
	Vector<dim, T> d[dim];

	FORCE_INLINE Matrix() {
		for (int i = 0; i < dim; ++i) {
			d[i] = std::move(Vector<dim, T>());
		}
	}

	FORCE_INLINE Matrix(const Matrix<dim, T> &m) {
		for (int i = 0; i < dim; ++i) {
			for (int j = 0; j < dim; ++j) {
				d[i][j] = m.d[i][j];
			}
		}
	}

	FORCE_INLINE Matrix(T a) {
		for (int i = 0; i < dim; ++i) {
			for (int j = 0; j < dim; ++j) {
				d[i][j] = a;
			}
		}
	}

	FORCE_INLINE Matrix(const Matrix &m) {
		*this = m;
	}

	FORCE_INLINE Matrix(const Vector<dim, T> &v) {
		for (int i = 0; i < dim; i++)
			this->d[i][i] = v[i];
	}

	template <
		typename F,
		typename = std::enable_if_t<std::is_convertible<
		F,
		std::function<Vector<dim__, T>(int)>>::value,
		int>>
		FORCE_INLINE explicit Matrix(const F &f) {
		for (int i = 0; i < dim; i++)
			this->d[i] = f(i);
	}

	template<typename = std::enable_if_t<dim == 2, int>>
	FORCE_INLINE Matrix(const Vector<dim, T> &v0, const Vector<dim, T> &v1) {
		d[0] = v0;
		d[1] = v1;
	}

	template<typename = std::enable_if_t<dim == 2, int>>
	FORCE_INLINE Matrix(T m00, T m01,
						T m10, T m11) {
		d[0][0] = m00; d[1][0] = m01;
		d[0][1] = m10; d[1][1] = m11;
	}

	template<typename = std::enable_if_t<dim == 3, int>>
	FORCE_INLINE Matrix(const Vector<dim, T> &v0, const Vector<dim, T> &v1, const Vector<dim, T> &v2) {
		d[0] = v0;
		d[1] = v1;
		d[2] = v2;
	}

	template<typename = std::enable_if_t<dim == 3, int>>
	FORCE_INLINE Matrix(T m00, T m01, T m02,
						T m10, T m11, T m12,
						T m20, T m21, T m22) {
		d[0][0] = m00; d[1][0] = m01; d[2][0] = m02;
		d[0][1] = m10; d[1][1] = m11; d[2][1] = m12;
		d[0][2] = m20; d[1][2] = m21; d[2][2] = m22;
	}

	template<typename = std::enable_if_t<dim == 4, int>>
	FORCE_INLINE Matrix(const Vector<dim, T> &v0, const Vector<dim, T> &v1, 
		const Vector<dim, T> &v2, const Vector<dim, T> &v3) {
		d[0] = v0;
		d[1] = v1;
		d[2] = v2;
		d[3] = v3;
	}

	template<typename = std::enable_if_t<dim == 3, int>>
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
			this->d[i] = std::move(m.d[i]);
		}
		return *this;
	}

	FORCE_INLINE Vector<dim, T>& operator[](int i) {
		return d[i];
	}

	FORCE_INLINE const T& operator()(int i, int j) {
		return d[j][i];
	}

	FORCE_INLINE T& operator()(int i, int j) {
		return d[j][i];
	}

	FORCE_INLINE Matrix operator+(const Matrix &a) const {
		return Matrix([=](int i) {this->d[i] + a.d[i]; });
	}

	FORCE_INLINE Matrix& operator+=(const Matrix &a) {
		for (int i = 0; i < dim; ++i) {
			this->d[i] += a.d[i];
		}
		return *this;
	}

	FORCE_INLINE Matrix operator-(const Matrix &a) const {
		return Matrix([=](int i) {this->d[i] - a.d[i]; });
	}

	FORCE_INLINE Matrix& operator-=(const Matrix &a) {
		for (int i = 0; i < dim; ++i) {
			this->d[i] -= a.d[i];
		}
		return *this;
	}

	FORCE_INLINE Vector<dim, T> operator*(const Vector<dim, T> &a) const {
		Vector<dim, T> ret(d[0] * a.d[0]);
		for (int i = 1; i < dim; ++i) {
			ret += d[i] * a.d[i];
		}
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
};




NAMESPACE_END