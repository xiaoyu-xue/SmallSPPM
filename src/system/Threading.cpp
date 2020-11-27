#include "Threading.h"

NAMESPACE_BEGIN

int ThreadIndex() {
	return tbb::task_arena::current_thread_index();
}

int NumSystemCores() {
	return std::max(1u, std::thread::hardware_concurrency());
}

int ThreadsNumber() {
	return NumSystemCores();
}

NAMESPACE_END
