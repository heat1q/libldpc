#include "cudamgd.h"

using namespace std;
using namespace ldpc;

cuda_mgd::cuda_mgd(const bool mgd) : mIsMgd(mgd) {}

void* cuda_mgd::operator new(size_t len)
{
	void* ptr;
	cudaMallocManaged(&ptr, len);
	cudaDeviceSynchronize();
	return ptr;
}

void cuda_mgd::operator delete(void* ptr)
{
	cudaDeviceSynchronize();
	cudaFree(ptr);
}
