#include "gtest/gtest.h"
#include "common/Platform.h"

int main(int argc, char* argv[]) {
	testing::InitGoogleTest(&argc, argv);
	RUN_ALL_TESTS();
}