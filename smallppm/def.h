#pragma once
#include <iostream>


#define ASSERT(expr) \
	do { if(!(expr)) { std::cerr << "Error: assertion `"#expr"' failed at " << __FILE__ << ":" << __LINE__ << std::endl; exit(2); } } while(0)
