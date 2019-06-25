#pragma once
#include <iostream>
#include <cassert>

/*
#define ASSERT(expr) \
	do { if(!(expr)) { std::cerr << "Error: assertion `"#expr"' failed at " << __FILE__ << ":" << __LINE__ << std::endl; exit(2); } } while(0)
*/

#define ASSERT assert
#define NAMESPACE_BEGIN namespace smallsppm {
#define NAMESPACE_END }