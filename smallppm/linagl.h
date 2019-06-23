#pragma once

#include "utils.h"

template<int dim, typename T>
struct Vector {

	Vector() {
		for (int i = 0; i < dim; ++i) {
			d[i] = T(0);
		}
	}

	Vector(const Vector<dim, T> &a) {
		for (int i = 0; i < dim; ++i) {
			d[i] = a.d[i];
		}
	}


	FORCE_INLINE Vector operator+(const Vector &b) const {
		Vector<dim, T> ret;
		for (int i = 0; i < dim; ++i) {
			ret.d[i] = d[i] + b.d[i];
		}
		return ret;
	}


	FORCE_INLINE Vector operator-(const Vector &b) const {
		Vector<dim, T> ret;
		for (int i = 0; i < dim; ++i) {
			ret.d[i] = d[i] - b.d[i];
		}
		return ret;
	}

	FORCE_INLINE Vector operator*(const Vector &b) const {
		Vector<dim, T> ret;
		for (int i = 0; i < dim; ++i) {
			ret.d[i] = d[i] * b.d[i];
		}
		return ret;
	}

	FORCE_INLINE Vector operator*(const T b) const {
		Vector<dim, T> ret;
		for (int i = 0; i < dim; ++i) {
			ret.d[i] = d[i] * b;
		}
		return ret;
	}

	FORCE_INLINE Vector operator/(const Vector &b) const {
		Vector<dim, T> ret;
		for (int i = 0; i < dim; ++i) {
			ret.d[i] = d[i] / b.d[i];
		}
		return ret;
	}

	FORCE_INLINE Vector operator/(const T b) const {
		Vector<dim, T> ret;
		for (int i = 0; i < dim; ++i) {
			ret.d[i] = d[i] / b;
		}
		return ret;
	}

	FORCE_INLINE bool operator==(const Vector &b) const {
		for (int i = 0; i < dim; ++i) {
			if (d[i] != b.[i]) {
				return false;
			}
		}
		return true;
	}

	FORCE_INLINE bool operator!=(const Vector &b) const {
		for (int i = 0; i < dim; ++i) {
			if (d[i] == b.[i]) {
				return false;
			}
		}
		return true;
	}

	FORCE_INLINE Vector& operator+=(const Vector &b) {
		for (int i = 0; i < dim; ++i) {
			d[i] += b.d[i];
		}
		return *this;
	}

	FORCE_INLINE Vector& operator-=(const Vector &b) {
		for (int i = 0; i < dim; ++i) {
			d[i] -= b.d[i];
		}
		return *this;
	}

	FORCE_INLINE Vector& operator*=(const Vector &b) {
		for (int i = 0; i < dim; ++i) {
			d[i] *= b.d[i];
		}
		return *this;
	}

	FORCE_INLINE Vector& operator/=(const Vector &b) {
		for (int i = 0; i < dim; ++i) {
			d[i] / = b.d[i];
		}
		return *this;
	}

	static constexpr int elements = dim;
	T d[dim];
};


template<typename T>
struct Vector<1, T> {
	static constexpr int elements = 1;
	union {
		T d[1];
		struct {
			T x;
		};
	};

	Vector(real x = 0){
		this->x = x;
	}

};

template<typename T>
struct Vector<2, T> {
	static constexpr int elements = 2;
	union {
		T d[2];
		struct {
			T x, y;
		};
	};

	Vector(real x = 0, real y = 0) {
		this->x = x;
		this->y = y;
	}
};

template<typename T>
struct Vector<3, T> {
	static constexpr int elements = 3;
	union {
		T d[3];
		struct {
			T x, y, z;
		};
	};

	Vector(real x = 0, real y = 0, real z = 0) {
		this->x = x;
		this->y = y;
		this->z = z;
	}
};

template<typename T>
struct Vector<4, T> {
	static constexpr int elements = 4;
	union {
		T d[4];
		struct {
			T x, y, z, w;
		};
	};

	Vector(real x = 0, real y = 0, real z = 0, real w = 0) {
		this->x = x;
		this->y = y;
		this->z = z;
		this->w = w;
	}
};


template<int dim, typename T>
FORCE_INLINE Vector<dim, T> operator*(const T a, const Vector<dim, T> &b) {
	Vector<dim, T> ret;
	for (int i = 0; i < dim; ++i) {
		ret.d[i] = b.d[i] * a;
	}
	return ret;
}

template<int dim, typename T>
FORCE_INLINE T Dot(const Vector<dim, T> &a, const Vector<dim, T> &b) {
	T ret = 0.0;
	for (int i = 0; i < dim; ++i) {
		ret += a.d[i] * b.d[i];
	}
	return ret;
}

template<typename T>
FORCE_INLINE Vector<3, T> Cross(const Vector<3, T> &a, const Vector<3, T> &b) {
	return Vector<3, T>(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
}