#pragma once
#include "linagl.h"
#include "ray.h"
#include "numeric_utils.h"
#include "AABB.h"

NAMESPACE_BEGIN

class Transform{
public:
	Transform() {
		mat = Matrix4::Identity();
		invMat = Matrix4::Identity();
	}

	Transform(
		real m00, real m01, real m02, real m03,
		real m10, real m11, real m12, real m13,
		real m20, real m21, real m22, real m23,
		real m30, real m31, real m32, real m33)
	{
		mat = Matrix4(
			m00, m01, m02, m03,
			m10, m11, m12, m13,
			m20, m21, m22, m23,
			m30, m31, m32, m33);

		invMat = Inverse(mat);
	}

	Transform(const real m[4][4]) {
		mat = Matrix4(
			m[0][0], m[0][1], m[0][2], m[0][3],
			m[1][0], m[1][1], m[1][2], m[1][3],
			m[2][0], m[2][1], m[2][2], m[2][3],
			m[3][0], m[3][1], m[3][2], m[3][3]);

		invMat = Inverse(mat);
	}

	Transform(const Matrix4 &m) {
		mat = m;
		invMat = Inverse(m);
	}

	Transform(const Matrix4 &m, const Matrix4 &invM) {
		mat = m;
		invMat = Inverse(invM);
	}

	friend Transform Inverse(const Transform &t) {
		return Transform(t.invMat, t.mat);
	}

	friend  Transform Transpose(const Transform &t) {
		return Transform(t.mat.Transpose(), t.invMat.Transpose());
	}

	const Matrix4 &GetMatrix() const {
		return mat;
	}

	const Matrix4 &GetInverseMatrix() const {
		return invMat;
	}

	FORCE_INLINE Vector3 operator()(const Vector3 &p) const;
	FORCE_INLINE Vector3 operator()(const Vector3 &p, Vector3 *pError) const;
	FORCE_INLINE Vector3 operator()(const Vector3 &p, const Vector3 &pError, 
		Vector3 *absError) const;
	FORCE_INLINE Vector3 TransformVector(const Vector3 &v) const;
	FORCE_INLINE Vector3 TransformVector(const Vector3 &v, Vector3 *absError) const;
	FORCE_INLINE Vector3 TransformVector(const Vector3 &v, const Vector3 &vError,
		Vector3 *absError) const;
	FORCE_INLINE Vector3 TransformNormal(const Vector3 &n) const;
	FORCE_INLINE Ray operator()(const Ray &r, Vector3 *oError,
		Vector3 *dError) const;
	FORCE_INLINE Ray operator()(const Ray &r) const;
	FORCE_INLINE Transform operator*(const Transform &t) const;
	FORCE_INLINE AABB operator()(const AABB &bound) const;
private:
	Matrix4 mat, invMat;
};

FORCE_INLINE Vector3 Transform::operator()(const Vector3 &p) const {
	Vector4 pp(p.x, p.y, p.z, 1.0);
	pp = mat * pp;
	if(pp.w == 1) {
		return Vector3(pp.x, pp.y, pp.z);
	}
	else {
		return Vector3(pp.x, pp.y, pp.z) / pp.w;
	}
}

FORCE_INLINE Vector3 Transform::operator()(const Vector3 &p, Vector3 *pError) const {
	Vector4 pp(p.x, p.y, p.z, 1.0);
	pp = mat * pp;
	real xAbsSum = (std::abs(mat(0, 0) * p.x) + std::abs(mat(0, 1) * p.y) +
		std::abs(mat(0, 2) * p.z) + std::abs(mat(0, 3)));
	real yAbsSum = (std::abs(mat(1, 0) * p.x) + std::abs(mat(1, 1) * p.y) +
		std::abs(mat(1, 2) * p.z) + std::abs(mat(1, 3)));
	real zAbsSum = (std::abs(mat(2, 0) * p.x) + std::abs(mat(2, 1) * p.y) +
		std::abs(mat(2, 2) * p.z) + std::abs(mat(2, 3)));
	*pError = gamma(3) * Vector3(xAbsSum, yAbsSum, zAbsSum);
	if(pp.w == 1) {
		return Vector3(pp.x, pp.y, pp.z);
	}
	else {
		return Vector3(pp.x, pp.y, pp.z) / pp.w;
	}
}

FORCE_INLINE Vector3 Transform::operator()(const Vector3 &p, const Vector3 &pError,
	Vector3 *absError) const {
	real x = p.x, y = p.y, z = p.z;
	Vector4 pp(p.x, p.y, p.z, 1.0);
	pp = mat * pp;
	absError->x =
		(gamma(3) + (real)1) *
		(std::abs(mat(0, 0)) * pError.x + std::abs(mat(0, 1)) * pError.y +
			std::abs(mat(0, 2)) * pError.z) +
		gamma(3) * (std::abs(mat(0, 0) * x) + std::abs(mat(0, 1) * y) +
			std::abs(mat(0, 2) * z) + std::abs(mat(0, 3)));
	absError->y =
		(gamma(3) + (real)1) *
		(std::abs(mat(1, 0)) * pError.x + std::abs(mat(1, 1)) * pError.y +
			std::abs(mat(1, 2)) * pError.z) +
		gamma(3) * (std::abs(mat(1, 0) * x) + std::abs(mat(1, 1) * y) +
			std::abs(mat(1, 2) * z) + std::abs(mat(1, 3)));
	absError->z =
		(gamma(3) + (real)1) *
		(std::abs(mat(2, 0)) * pError.x + std::abs(mat(2, 1)) * pError.y +
			std::abs(mat(2, 2)) * pError.z) +
		gamma(3) * (std::abs(mat(2, 0) * x) + std::abs(mat(2, 1) * y) +
			std::abs(mat(2, 2) * z) + std::abs(mat(2, 3)));

	if (pp.w == 1.) {
		return Vector3(pp.x, pp.y, pp.z);
	}
	else {
		return Vector3(pp.x, pp.y, pp.z) / pp.w;
	}

}


FORCE_INLINE Vector3 Transform::TransformVector(const Vector3 &v) const {
	Vector4 vv(v.x, v.y, v.z, 0);
	vv = mat * vv;
	return Vector3(vv.x, vv.y, vv.z);
}

FORCE_INLINE Vector3 Transform::TransformVector(const Vector3 &v, Vector3 *absError) const {
	real x = v.x, y = v.y, z = v.z;
	absError->x =
		gamma(3) * (std::abs(mat(0, 0) * v.x) + std::abs(mat(0, 1) * v.y) +
			std::abs(mat(0, 2) * v.z));
	absError->y =
		gamma(3) * (std::abs(mat(1, 0) * v.x) + std::abs(mat(1, 1) * v.y) +
			std::abs(mat(1, 2) * v.z));
	absError->z =
		gamma(3) * (std::abs(mat(2, 0) * v.x) + std::abs(mat(2, 1) * v.y) +
			std::abs(mat(2, 2) * v.z));
	Vector4 vv(v.x, v.y, v.z, 0);
	vv = mat * vv;
	return Vector3(vv.x, vv.y, vv.z);
}

FORCE_INLINE Vector3 Transform::TransformVector(const Vector3 &v, const Vector3 &vError,
	Vector3 *absError) const {
	real x = v.x, y = v.y, z = v.z;
	absError->x =
		(gamma(3) + (real)1) *
		(std::abs(mat(0, 0)) * vError.x + std::abs(mat(0, 1)) * vError.y +
			std::abs(mat(0, 2)) * vError.z) +
		gamma(3) * (std::abs(mat(0, 0) * v.x) + std::abs(mat(0, 1) * v.y) +
			std::abs(mat(0, 2) * v.z));
	absError->y =
		(gamma(3) + (real)1) *
		(std::abs(mat(1, 0)) * vError.x + std::abs(mat(1, 1)) * vError.y +
			std::abs(mat(1, 2)) * vError.z) +
		gamma(3) * (std::abs(mat(1, 0) * v.x) + std::abs(mat(1, 1) * v.y) +
			std::abs(mat(1, 2) * v.z));
	absError->z =
		(gamma(3) + (real)1) *
		(std::abs(mat(2, 0)) * vError.x + std::abs(mat(2, 1)) * vError.y +
			std::abs(mat(2, 2)) * vError.z) +
		gamma(3) * (std::abs(mat(2, 0) * v.x) + std::abs(mat(2, 1) * v.y) +
			std::abs(mat(2, 2) * v.z));
	Vector4 vv(v.x, v.y, v.z, 0);
	vv = mat * vv;
	return Vector3(vv.x, vv.y, vv.z);
}

FORCE_INLINE Ray Transform::operator()(const Ray& r) const {
	Vector3 oError;
	Vector3 o = (*this)(r.o, &oError);
	Vector3 d = (*this)(r.d);

	real length2 = d.Length2();
	real tMax = r.tMax, tMin = r.tMin;
	if(length2 > 0) {
		real dt = Dot(Abs(d), oError) / length2;
		o += d * dt;
		tMax -= dt;
		tMin += dt;
	}
	return Ray(o, d, tMax, tMin);
}

FORCE_INLINE Ray Transform::operator()(const Ray &r, Vector3 *oError, Vector3 *dError) const {
	Vector3 o = (*this)(r.o, oError);
	Vector3 d = (*this)(r.d, dError);
	real tMax = r.tMax;
	real tMin = r.tMin;
	real length2 = d.Length2();
	if (length2 > 0) {
		real dt = Dot(Abs(d), *oError) / length2;
		o += d * dt;
		//        tMax -= dt;
	}
	return Ray(o, d, tMax, tMin);
}

FORCE_INLINE Vector3 Transform::TransformNormal(const Vector3 &n) const {
	Vector4 nn(n.x, n.y, n.z, 0);
	nn = invMat.Transpose() * nn;
	return Vector3(nn.x, nn.y, nn.z);
}

FORCE_INLINE Transform Transform::operator*(const Transform &t) const {
	return Transform(mat * t.mat, t.invMat * invMat);
}

FORCE_INLINE AABB Transform::operator()(const AABB &bound) const {
	return AABB((*this)(bound.minPoint), (*this)(bound.maxPoint));
}

NAMESPACE_END