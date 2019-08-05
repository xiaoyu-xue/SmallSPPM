#include "linagl.h"
#include "gtest/gtest.h"
#include "math_utils.h"


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
			mat1(1.1, 2, 3, 4,
				5, 6, 7.3, 8,
				9, 10, 11, 12,
				13, 14.9, 15, 16);
		Matrix<4, real, IntrinsicSet::None>
			mat2(1.1, 2, 3, 4,
				5, 6, 7.3, 8,
				9, 10, 11, 12,
				13, 14.9, 15, 16);
		EXPECT_EQ(Equal(mat1, mat2), true);
	}
};