#include "LowDiscrepency.h"

GY_NAMESPACE_BEGIN

int IsPrime(int a) noexcept {
	ASSERT(a >= 2);
	for (int i = 2; i * i <= a; i++) {
		if (a % i == 0)
			return false;
	}
	return true;
}

GY_NAMESPACE_END