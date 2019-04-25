#pragma once

#include <iostream>
#include <stdint.h>
#include <exception>

#define NUM_THREADS 128
#define gpuErrchk(ans) { ldpc::gpuAssert((ans), __FILE__, __LINE__); }

namespace ldpc
{
	class cuda_mgd
	{
	public:
		cuda_mgd(const bool mgd);

		void* operator new(size_t len);
		void operator delete(void* ptr);

	protected:
		bool mIsMgd = false;
	};

	inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
	{
		if (code != cudaSuccess)
		{
			fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
			if (abort) exit(code);
		}
	}
}
