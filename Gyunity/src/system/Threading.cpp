#include "Threading.h"

GYT_NAMESPACE_BEGIN

int ThreadIndex() {
	return tbb::task_arena::current_thread_index();
}

int NumSystemCores() {
	return std::max(1u, std::thread::hardware_concurrency());
}

int ThreadsNumber() {
	return NumSystemCores();
}

GYT_NAMESPACE_END
