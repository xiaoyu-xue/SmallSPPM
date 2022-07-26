#include "gtest/gtest.h"
#include "math/Linagl.h"
#include "math/MathUtils.h"

GYT_NAMESPACE_BEGIN
TEST(TestMatrixIsEqual, IsEqual) {
	{
		Matrix<4, real, IntrinsicSet::None>
			mat1(1, 2, 3, 4,
				5, 6, 7, 8,
				9, 10, 11, 12,
				13, 14, 15, 16);
		Matrix<4, real, IntrinsicSet::None>
			mat2(1, 2, 3, 4,
				5, 6, 7, 8,
				9, 10, 11, 12,
				13, 14, 15, 16);
		EXPECT_EQ(mat1 == mat2, true);
	}
	{
		Matrix<4, real, IntrinsicSet::None>
			mat1(1.1f, 2, 3, 4,
				5, 6, 7.3f, 8,
				9, 10, 11, 12,
				13, 14.9f, 15, 16);
		Matrix<4, real, IntrinsicSet::None>
			mat2(1.1f, 2, 3, 4,
				5, 6, 7.3f, 8,
				9, 10, 11, 12,
				13, 14.9f, 15, 16);
		EXPECT_EQ(Equal(mat1, mat2), true);
	}
};

TEST(TestMatrixNotEqual, NotEqual) {
	{
		Matrix<4, real, IntrinsicSet::None>
			mat1(1.1f, 2, 3, 4,
				5, 6, 7.3f, 8,
				9, 10, 11, 12,
				13, 14.9f, 15, 16);
		Matrix<4, real, IntrinsicSet::None>
			mat2(1.1f, 2, 3, 4,
				5, 6, 7.3f, 8,
				9, 10, 11, 13,
				13, 14.9f, 15, 16);
		EXPECT_EQ(mat1 != mat2, true);
	}
};

TEST(TestMatrixPlus, Plus) {
	{
		Matrix<4, real, IntrinsicSet::None>
			mat1(1, 2, 3, 4,
				5, 6, 7, 8,
				9, 10, 11, 12,
				13, 14, 15, 16);
		Matrix<4, real, IntrinsicSet::None>
			mat2(1, 2, 3, 4,
				5, 6, 7, 8,
				9, 10, 11, 12,
				13, 14, 15, 16);
		Matrix<4, real, IntrinsicSet::None>
			mat3(2, 4, 6, 8,
				10, 12, 14, 16,
				18, 20, 22, 24,
				26, 28, 30, 32);
		EXPECT_EQ(mat1 + mat2 == mat3, true);
	}
	{
		Matrix<4, real, IntrinsicSet::None>
			mat1(1, 2, 3, 4,
				5, 6, 7, 8,
				9, 10, 11, 12,
				13, 14, 15, 16);
		Matrix<4, real, IntrinsicSet::None>
			mat2(1, 2, 3, 4,
				5, 6, 7, 8,
				9, 10, 11, 12,
				13, 14, 15, 16);
		Matrix<4, real, IntrinsicSet::None>
			mat3(2, 4, 6, 8,
				10, 12, 14, 16,
				18, 20, 22, 24,
				26, 28, 30, 32);
		mat1 += mat2;
		EXPECT_EQ(mat1 == mat3, true);
	}
};

TEST(TestMatrixMinus, Minus) {
	{
		Matrix<4, real, IntrinsicSet::None>
			mat1(1, 2, 3, 4,
				5, 6, 7, 8,
				9, 10, 11, 12,
				13, 14, 15, 16);
		Matrix<4, real, IntrinsicSet::None>
			mat2(1, 2, 3, 4,
				5, 6, 7, 8,
				9, 10, 11, 12,
				13, 14, 15, 16);
		Matrix<4, real, IntrinsicSet::None>
			mat3 = Matrix<4, real, IntrinsicSet::None>::Zeros();
		EXPECT_EQ((mat1 - mat2) == mat3, true);
	}
	{
		Matrix<4, real, IntrinsicSet::None>
			mat1(1, 2, 3, 4,
				5, 6, 7, 8,
				9, 10, 11, 12,
				13, 14, 15, 16);
		Matrix<4, real, IntrinsicSet::None>
			mat2(1, 2, 3, 4,
				5, 6, 7, 8,
				9, 10, 11, 12,
				13, 14, 15, 16);
		Matrix<4, real, IntrinsicSet::None>
			mat3 = Matrix<4, real, IntrinsicSet::None>::Zeros();
		mat1 -= mat2;
		EXPECT_EQ(mat1 == mat3, true);
	}
};


TEST(TestMatrixMulMatrix, MulMatrix) {
	{
		Matrix<4, real, IntrinsicSet::None>
			mat1(1, 2, 3, 4,
				5, 6, 7, 8,
				9, 10, 11, 12,
				13, 14, 15, 16);
		Matrix<4, real, IntrinsicSet::None>
			mat2(1, 2, 3, 4,
				5, 6, 7, 8,
				3, 2, 1, 6,
				1, 3, 5, 7);
		Matrix<4, real, IntrinsicSet::None>
			mat3(24, 32, 40, 66,
				 64, 84, 104, 166,
				 104, 136, 168, 266,
				 144, 188, 232, 366);
		EXPECT_EQ(mat1 * mat2 == mat3, true);
	}
	{
		Matrix<4, real, IntrinsicSet::None>
			mat1(1, 2, 3, 4,
				5, 6, 7, 8,
				9, 10, 11, 12,
				13, 14, 15, 16);
		Matrix<4, real, IntrinsicSet::None>
			mat2(1, 2, 3, 4,
				5, 6, 7, 8,
				3, 2, 1, 6,
				1, 3, 5, 7);
		Matrix<4, real, IntrinsicSet::None>
			mat3(24, 32, 40, 66,
				64, 84, 104, 166,
				104, 136, 168, 266,
				144, 188, 232, 366);
		mat1 *= mat2;
		EXPECT_EQ(mat1 == mat3, true);
	}
	{
		Matrix<4, real, IntrinsicSet::None>
			mat1(1, 2, 3, 4,
				5, 6, 7, 8,
				9, 10, 11, 12,
				13, 14, 15, 16);
		Matrix<4, real, IntrinsicSet::None>
			mat2(2, 4, 6, 8,
				10, 12, 14, 16,
				18, 20, 22, 24,
				26, 28, 30, 32);

		EXPECT_EQ(mat1 * 2 == mat2, true);
		EXPECT_EQ(2 * mat1 == mat2, true);
	}
	{
		Matrix<4, real, IntrinsicSet::None>
			mat1(1, 2, 3, 4,
				5, 6, 7, 8,
				9, 10, 11, 12,
				13, 14, 15, 16);
		Matrix<4, real, IntrinsicSet::None>
			mat2(2, 4, 6, 8,
				10, 12, 14, 16,
				18, 20, 22, 24,
				26, 28, 30, 32);
		mat1 *= 2;
		EXPECT_EQ(mat1 == mat2, true);

	}
};


TEST(TestMatrixMulVector, MulVector) {
	{
		Matrix<4, real, IntrinsicSet::None>
			mat(1, 2, 3, 4,
				5, 6, 7, 8,
				9, 10, 11, 12,
				13, 14, 15, 16);
		Vector<4, real, IntrinsicSet::None> v1(1.1f, 2.2f, 3.3f, 4.4f);
		Vector<4, real, IntrinsicSet::None> v2(33.f, 77.f, 121.f, 165.f);
		EXPECT_EQ(Equal(mat * v1, v2), true);
	}
	{
		Matrix<4, real, IntrinsicSet::None>
			mat(0.992797f, 0.307158f, 0.878931f, 0.862293f,
				0.688249f, 0.42504f, 0.677657f, 0.540247f,
				0.970518f, 0.0388582f, 0.482551f, 0.547493f,
				0.973952f, 0.116856f, 0.645318f, 0.350542f);
		Vector<4, real, IntrinsicSet::None> 
		v1(0.8310123813548647f,
			0.7829980173878517f,
			0.48269930108431525f,
			0.4714602712916691f);
		Vector<4, real, IntrinsicSet::None> 
		v2(1.8963264168581206f,
			1.4865590956869474f,
			1.3279869317662532f,
			1.3776250350210677f);
		EXPECT_EQ(Equal(mat * v1, v2), true);
	}
};


TEST(TestMatrixDeterminant, Determinant) {
	{
		Matrix<4, real, IntrinsicSet::None>
			mat(1, 2, 3, 4,
				5, 6, 7, 8,
				9, 10, 11, 12,
				13, 14, 15, 16);
		EXPECT_EQ(Determinant(mat) == 0, true);
	}
	{
		Matrix<4, real, IntrinsicSet::None>
			mat(0.992797f, 0.307158f, 0.878931f, 0.862293f,
				0.688249f, 0.42504f, 0.677657f, 0.540247f,
				0.970518f, 0.0388582f, 0.482551f, 0.547493f,
				0.973952f, 0.116856f, 0.645318f, 0.350542f);
		real expectedDeterminant = 0.02143648345371967f;
		EXPECT_EQ(Equal(Determinant(mat), expectedDeterminant), true);
	}
};

TEST(TestMatrixInv, Inv) {
	{
		Matrix<4, real, IntrinsicSet::None>
			scaleMatrix(1, 0, 0, 0,
				0, 2, 0, 0,
				0, 0, 5, 0,
				0, 0, 0, 1);
		Matrix<4, real, IntrinsicSet::None>
			invScaleMatrix(1.0f, 0, 0, 0,
						   0, 0.5f, 0, 0,
						   0, 0, 0.2f, 0,
						   0, 0, 0, 1.0f);
		EXPECT_EQ(Equal(Inverse(scaleMatrix), invScaleMatrix), true);
	}
	{
		Matrix<4, real, IntrinsicSet::None>
			translateMatrix(1, 0, 0, 1.1f,
							0, 1, 0, 2.2f,
							0, 0, 1, 3.3f,
							0, 0, 0, 1);
		Matrix<4, real, IntrinsicSet::None>
			invTranslateMatrix(1.0f, 0, 0, -1.1f,
							   0, 1.0f, 0, -2.2f,
							   0, 0, 1.0f, -3.3f,
							   0, 0, 0, 1.0f);
		EXPECT_EQ(Equal(Inverse(translateMatrix), invTranslateMatrix), true);
	}
	{
		Matrix<4, real, IntrinsicSet::None>
			mat(0.992797f, 0.307158f, 0.878931f, 0.862293f,
				0.688249f, 0.42504f, 0.677657f, 0.540247f,
				0.970518f, 0.0388582f, 0.482551f, 0.547493f,
				0.973952f, 0.116856f, 0.645318f, 0.350542f);
		Matrix<4, real, IntrinsicSet::None>
			expectedInversedMat(-2.84865f, 1.83357f, 2.73568f, -0.091224f,
						-4.12886f, 5.67329f, 2.146f, -1.93875f,
						4.11491f, -3.47325f, -5.35782f, 3.5988f,
						1.71591f, -0.591686f, 1.54706f, -2.8726f);
		EXPECT_EQ(Equal(Inverse(mat), expectedInversedMat, 1e-4), true);
	}
};

GYT_NAMESPACE_END