#include <gtest/gtest.h>
#include <iostream>

int main(int argc, char* argv[])
{
	testing::InitGoogleTest(&argc, argv);
	//    testing::GTEST_FLAG(filter) = "SharedPtr.Struct";
	int r = RUN_ALL_TESTS();
	std::cin.get();
	return r;
}
