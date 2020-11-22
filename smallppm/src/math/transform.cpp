#include "transform.h"
#include "math/math_utils.h"

Transform Transform::Translate(const Vector3 &v) {
	Matrix4 m(1, 0, 0, v.x, 
			  0, 1, 0, v.y, 
			  0, 0, 1, v.z, 
		      0, 0, 0, 1);
	Matrix4 minv(1, 0, 0, -v.x, 
				 0, 1, 0, -v.y, 
				 0, 0, 1, -v.z, 
				 0, 0, 0, 1);
	return Transform(m, minv);
}

Transform Transform::Scale(real sx, real sy, real sz) {
	Matrix4 m(sx, 0, 0, 0, 
			  0, sy, 0, 0, 
			  0, 0, sz, 0, 
			  0, 0, 0, 1);
	Matrix4 minv(1 / sx, 0, 0, 0, 
				 0, 1 / sy, 0, 0, 
				 0, 0, 1 / sz, 0, 
				 0, 0, 0, 1);
	return Transform(m, minv);
}

Transform Transform::RotateX(real theta) {
	real sinTheta = std::sin(Radians(theta));
	real cosTheta = std::cos(Radians(theta));
	Matrix4 m(1, 0, 0, 0, 
			  0, cosTheta, -sinTheta, 0, 
			  0, sinTheta, cosTheta, 0,
			  0, 0, 0, 1);
	return Transform(m, m.Transpose());
}

Transform Transform::RotateY(real theta) {
	real sinTheta = std::sin(Radians(theta));
	real cosTheta = std::cos(Radians(theta));
	Matrix4 m(cosTheta, 0, sinTheta, 0, 
			  0, 1, 0, 0, 
			  -sinTheta, 0, cosTheta, 0,
			  0, 0, 0, 1);
	return Transform(m, m.Transpose());
}

Transform Transform::RotateZ(real theta) {
	real sinTheta = std::sin(Radians(theta));
	real cosTheta = std::cos(Radians(theta));
	Matrix4 m(cosTheta, -sinTheta, 0, 0, 
			  sinTheta, cosTheta, 0, 0, 
			  0, 0, 1, 0,
			  0, 0, 0, 1);
	return Transform(m, m.Transpose());
}

Transform Transform::Rotate(real theta, const Vector3 &axis) {
	Vector3 a = axis.Norm();
	real sinTheta = std::sin(Radians(theta));
	real cosTheta = std::cos(Radians(theta));

	Matrix4 mat(a.x * a.x + (1 - a.x * a.x) * cosTheta,
				a.x * a.y * (1 - cosTheta) - a.z * sinTheta,
				a.x * a.z * (1 - cosTheta) + a.y * sinTheta,
				0,
				a.x * a.y * (1 - cosTheta) + a.z * sinTheta,
				a.y * a.y + (1 - a.y * a.y) * cosTheta,
				a.y * a.z * (1 - cosTheta) - a.x * sinTheta,
				0,
				a.x * a.z * (1 - cosTheta) - a.y * sinTheta,
				a.y * a.z * (1 - cosTheta) + a.x * sinTheta,
				a.z * a.z + (1 - a.z * a.z) * cosTheta,
				0,
				0, 0, 0, 1);
	return Transform(mat, mat.Transpose());
}

Transform Transform::LookAt(const Vector3 &pos, const Vector3 &look, const Vector3 &up) {
	Matrix4 cameraToWorld;
	// Initialize fourth column of viewing matrix
	cameraToWorld(0, 3) = pos.x;
	cameraToWorld(1, 3) = pos.y;
	cameraToWorld(2, 3) = pos.z;
	cameraToWorld(3, 3) = 1;

	// Initialize first three columns of viewing matrix
	Vector3 dir = (look - pos).Norm();
	//if (Cross(up.Norm(), dir).Length() == 0) {
	//	Error(
	//		"\"up\" vector (%f, %f, %f) and viewing direction (%f, %f, %f) "
	//		"passed to LookAt are pointing in the same direction.  Using "
	//		"the identity transformation.",
	//		up.x, up.y, up.z, dir.x, dir.y, dir.z);
	//	return Transform();
	//}


	Vector3 right = Cross(up.Norm(), dir).Norm();
	Vector3 newUp = Cross(dir, right);
	cameraToWorld(0, 0) = right.x;
	cameraToWorld(1, 0) = right.y;
	cameraToWorld(2, 0) = right.z;
	cameraToWorld(3, 0) = 0.;
	cameraToWorld(0, 1) = newUp.x;
	cameraToWorld(1, 1) = newUp.y;
	cameraToWorld(2, 1) = newUp.z;
	cameraToWorld(3, 1) = 0.;
	cameraToWorld(0, 2) = dir.x;
	cameraToWorld(1, 2) = dir.y;
	cameraToWorld(2, 2) = dir.z;
	cameraToWorld(3, 2) = 0.;
	return Transform(Inverse(cameraToWorld), cameraToWorld);
}

Transform Transform::Orthographic(real zNear, real zFar) {
	return Scale(1, 1, 1 / (zFar - zNear)) * Translate(Vector3(0, 0, -zNear));
}

Transform Transform::Perspective(real fovy, real aspect, real dis, real n, real f) {

	real t = dis * std::tan(Radians(fovy) * 0.5);
	real r = t * aspect;
	//NDC x in [-1, 1], y in [-1, 1], z in [-1, 1]
	//Matrix4 persp(dis / r, 0, 0, 0,
	//			  0, dis / t, 0, 0,
	//			  0, 0, (f + n) / (f - n), -2 * f * n / (f - n),
	//			  0, 0, 1, 0);

	//NDC x in [-1, 1], y in [-1, 1], z in [0, 1]
	//We use unstandard NDC since z in [0, 1] is more natural
	Matrix4 persp(dis / r, 0, 0, 0,
				0, dis / t, 0, 0,
				0, 0, f / (f - n), - f * n / (f - n),
				0, 0, 1, 0);
	return Transform(persp, Inverse(persp));
}