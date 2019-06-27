#include "linagl.h"
#include "gtest/gtest.h"
#include "math_utils.h"

TEST(TestOstreamOperator, Ostream) {
	{
		std::ostringstream os;
		os << Vector3(1, 2, 3);
		EXPECT_EQ(os.str(), "[ 1, 2, 3 ]");
	}
	{
		std::ostringstream os;
		os << Vector3(1.3, 2.5, 3.7);
		EXPECT_EQ(os.str(), "[ 1.3, 2.5, 3.7 ]");
	}
};

TEST(TestOperatorIsEqual, IsEqual) {
	{
		EXPECT_EQ(Vector3(1.1, 2.2, 3.3) == Vector3(1.1, 2.2, 3.3), true);
	}
	{
		EXPECT_EQ(Vector3(1.1, 2.2, 3.5) == Vector3(1.1, 2.2, 3.3), false);
	}
}

TEST(TestOperatorEqual, Equal) {
	{
		Vector3 a, b(1, 2, 3);
		a = b;
		EXPECT_EQ(a == b, true);
	}
}

TEST(TestOperatorNotEqual, NotEqual) {
	{
		EXPECT_EQ(Vector3(1.1, 2.2, 3.3) != Vector3(1.1, 2.2, 3.3), false);
	}
	{
		EXPECT_EQ(Vector3(1.1, 2.2, 3.5) != Vector3(1.1, 2.2, 3.3), true);
	}
}

TEST(TestOperatorPlus, Plus) {
	{
		Vector3 a(1, 2, 3), b(3, 4, 5), c(4, 6, 8);
		EXPECT_EQ(a + b == c, true);
	}
	{
		Vector3 a(1, 2, 3), b(3, 4, 5), c(4, 6, 8);
		a += b;
		EXPECT_EQ(a == c, true);
	}
}

TEST(TestOperatorMinus, Minus) {
	{
		Vector3 a(1, 2, 3), b(3, 4, 5), c(-2, -2, -2);
		EXPECT_EQ(a - b == c, true);
	}
	{
		Vector3 a(1, 2, 3), b(3, 4, 5), c(-2, -2, -2);
		a -= b;
		EXPECT_EQ(a == c, true);
	}
	{
		Vector3 a(1, 2, 3), b(-1, -2, -3);
		EXPECT_EQ(-a == b, true);
	}
}

TEST(TestOperatorMul, Mul) {
	{
		Vector3 a(1, 2, 3), b(3, 4, 5), c(3, 8, 15);
		EXPECT_EQ(a * b == c, true);
	}
	{
		Vector3 a(1, 2, 3), c(3, 6, 9);
		real b = 3;
		EXPECT_EQ(a * b == c, true);
	}
	{
		Vector3 a(1, 2, 3), c(3, 6, 9);
		real b = 3;
		EXPECT_EQ(b * a == c, true);
	}
	{
		Vector3 a(1, 2, 3), b(3, 4, 5), c(3, 8, 15);
		a *= b;
		EXPECT_EQ(a == c, true);
	}
}


TEST(TestOperatorDiv, Div) {
	{
		Vector3 a(3, 8, 15), b(1, 2, 3), c(3, 4, 5);
		EXPECT_EQ(a / b == c, true);
	}
	{
		Vector3 a(3, 6, 9), c(1, 2, 3);
		real b = 3;
		EXPECT_EQ(a / b == c, true);
	}
	{
		Vector3 a(3, 8, 15), b(1, 2, 3), c(3, 4, 5);
		a /= b;
		EXPECT_EQ(a == c, true);
	}
}

TEST(TestDot, Dot) {
	{
		Vector3 a(1, 2, 3), b(5, 6, 7);
		real c = 38;
		EXPECT_EQ(a.Dot(b), c);
	}
	{
		Vector3 a(1, 2, 3), b(5, 6, 7);
		real c = 38;
		EXPECT_EQ(b.Dot(a), c);
	}
	{
		Vector3 a(1, 2, 3), b(5, 6, 7);
		real c = 38;
		EXPECT_EQ(Dot(a, b), c);
	}
	{
		Vector3 a(1, 2, 3), b(5, 6, 7);
		real c = 38;
		EXPECT_EQ(Dot(b, a), c);
	}
}

TEST(TestCross, Cross) {
	{
		Vector3 a(1, 2, 3), b(4, 5, 6), c(-3, 6, -3);
		EXPECT_EQ(a.Cross(b) == c, true);
	}
	{
		Vector3 a(1, 2, 3), b(4, 5, 6), c(-3, 6, -3);
		EXPECT_EQ(Cross(a, b) == c, true);
	}
	{
		Vector3 a(1, 2, 3), b(4, 5, 6), c(3, -6, 3);
		EXPECT_EQ(b.Cross(a) == c, true);
	}
	{
		Vector3 a(1, 2, 3), b(4, 5, 6), c(3, -6, 3);
		EXPECT_EQ(Cross(b, a) == c, true);
	}
	{
		Vector3 a(1, 2, 3), b(4, 5, 6), c(3, -6, 3);
		EXPECT_EQ(a.Cross(b) == -b.Cross(a), true);
	}
	{
		Vector3 a(1, 2, 3), b(4, 5, 6), c(3, -6, 3);
		EXPECT_EQ(Cross(a, b) == -Cross(b, a), true);
	}
}

TEST(TestLength, Length) {
	{
		Vector3 a(1, 2, 3);
		real length2 = 14;
		EXPECT_EQ(a.Length2(), length2);
	}
	{
		Vector3 a(1, 2, 3);
		real length = 3.741657386;
		EXPECT_EQ(Equal(a.Length(), length), true);
	}
}

TEST(TestNormalize, Normalize) {
	{
		Vector3 a(1, 2, 3), b(0.2672612419, 0.5345224838, 0.80178372573);
		EXPECT_EQ(Equal(a.Normal(), b), true);
	}
	{
		Vector3 a(1, 2, 3), b(0.2672612419, 0.5345224838, 0.80178372573);
		a.Normalize();
		EXPECT_EQ(Equal(a, b), true);
	}
	{
		Vector3 a(1, 2, 3), b(0.2672612419, 0.5345224838, 0.80178372573);
		Normalize(a);
		EXPECT_EQ(Equal(a, b), true);
	}
}