#include "linagl.h"
#include "gtest/gtest.h"
#include "math_utils.h"
#include "transform.h"
#include "pch.h"


TEST(TransformTest, Translate) {
	{
		Vec3 vector3(1, 2, 3);
		Vec3 point3(1, 2, 3);
		Vec3 normal(1, 2, 3);
		Ray ray(Vec3(1, 2, 3), Vec3(1, 0, 0));
		Transform transform = Transform::Translate(Vec3(1, 2, 3));

		Vec3 absError;
		Vec3 transformedVector = transform.TransformVector(vector3, &absError);
		Vec3 transformedPoint = transform(point3, &absError);
		EXPECT_EQ(transformedVector, Vec3(4, 5, 6));
		//EXPECT_EQ()
	}
}