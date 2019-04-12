#include "../ldpcsim.h"

#include <math.h>
#include <exception>

using namespace ldpc;
using namespace thrust;


constellation::constellation(const uint16_t pM)
	: mM(pM), mLog2M(log2(pM)), cuda_mgd(true)
{
	mPX = device_vector<double>(mM, 1.0/mM);
	mX = device_vector<double>(mM);

	double m = 0;
	for (size_t j = 0; j < mM; ++j)
	{
		mX[j] = (double) -mM+1+2*j;
		m += mX[j] * mX[j] * mPX[j];
	}

	for (size_t j = 0; j < mM; ++j) {
		mX[j] = mX[j]/sqrt(m);
	}
}
