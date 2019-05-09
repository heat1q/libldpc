#include "ldpc.h"


namespace ldpc
{
/*
* Code device
*/
//init constructor
ldpc_code::ldpc_code(const char *pFileName)
	: mMaxDC(0)
{
	try
	{
		FILE *fpCode = fopen(pFileName, "r");
		if (!fpCode)
		{
			throw std::runtime_error("can not open codefile for reading.");
		}

		fscanf(fpCode, "nc: %lu\n", &mN);
		fscanf(fpCode, "mc: %lu\n", &mM);
		fscanf(fpCode, "nct: %lu\n", &mNCT);
		fscanf(fpCode, "mct: %lu\n", &mMCT);
		fscanf(fpCode, "nnz: %lu\n", &mNNZ);
		mK = mN - mM;
		mKCT = mNCT - mMCT;

		fscanf(fpCode, "puncture [%lu]: ", &(mNumPuncture));
		mNumPunctureSys = 0;
		mNumPuncturePar = 0;
		if (mNumPuncture != 0)
		{
			mPuncture = vec_size_t(mNumPuncture);
			for (std::size_t i = 0; i < mNumPuncture; i++)
			{
				fscanf(fpCode, " %lu ", &(mPuncture[i]));
				if (mPuncture[i] < mK)
				{
					mNumPunctureSys++;
				}
				else
				{
					mNumPuncturePar++;
				}
			}
		}

		fscanf(fpCode, "shorten [%lu]: ", &mNumShorten);
		if (mNumShorten != 0)
		{
			mShorten = vec_size_t(mNumShorten);
			for (std::size_t i = 0; i < mNumShorten; i++)
			{
				fscanf(fpCode, " %lu ", &(mShorten[i]));
			}
		}

		vec_size_t cwTmp(mM);
		vec_size_t vwTmp(mN);

		mCW = vec_size_t(mM);
		mVW = vec_size_t(mN);
		mR = vec_size_t(mNNZ);
		mC = vec_size_t(mNNZ);

		for (std::size_t i = 0; i < mNNZ; i++)
		{
			fscanf(fpCode, "%lu %lu\n", &(mR[i]), &(mC[i]));
			mCW[mR[i]]++;
			mVW[mC[i]]++;
		}

		mCN = mat_size_t(mM, vec_size_t());
		for (std::size_t i = 0; i < mM; i++)
		{
			mCN[i] = vec_size_t(mCW[i]);
		}

		mVN = mat_size_t(mN, vec_size_t());
		for (std::size_t i = 0; i < mN; i++)
		{
			mVN[i] = vec_size_t(mVW[i]);
		}

		for (std::size_t i = 0; i < mNNZ; i++)
		{
			mCN[mR[i]][cwTmp[mR[i]]++] = i;
			mVN[mC[i]][vwTmp[mC[i]]++] = i;
		}

		for (std::size_t i = 0; i < mM; i++)
		{
			if (mCW[i] > mMaxDC)
			{
				mMaxDC = mCW[i];
			}
		}

		fclose(fpCode);
	}
	catch (std::exception &e)
	{
		std::cout << "Error:   ldpc_code::  ldpc_code(): " << e.what() << std::endl;
		exit(EXIT_FAILURE);
	}
}

void ldpc_code::print()
{
	std::cout << "nc : " << mN << "\n";
	std::cout << "mc : " << mM << "\n";
	std::cout << "kc : " << mK << "\n";
	std::cout << "nnz : " << mNNZ << "\n";
	std::cout << "nct :" << mNCT << "\n";
	std::cout << "mct : " << mMCT << "\n";
	std::cout << "kct : " << mKCT << "\n";
	std::cout << "max dc : " << mMaxDC << "\n";
	std::cout << "num puncture: " << mNumPuncture << "\n";
	std::cout << "num puncture sys: " << mNumPunctureSys << "\n";
	std::cout << "num puncture par: " << mNumPuncturePar << "\n";
	std::cout << "num shorten: " << mNumShorten << "\n";
}

void dec2bin(std::size_t val, uint8_t m)
{
	for (std::size_t i = 0; i < m; i++)
		printf("%lu", (val >> (m - i - 1) & 0x01));
}

double jacobian(double L1, double L2)
{
#ifdef CN_APPROX_LIN
	return sign(L1) * sign(L2) * fmin(fabs(L1), fabs(L2)) + jacobian_lin_approx(L1 + L2) - jacobian_lin_approx(L1 - L2);
#elif defined CN_APPROX_MINSUM
	return sign(L1) * sign(L2) * fmin(fabs(L1), fabs(L2));
#else
	return sign(L1) * sign(L2) * fmin(fabs(L1), fabs(L2)) + log((1 + exp(-fabs(L1 + L2))) / (1 + exp(-fabs(L1 - L2))));
#endif
}

double jacobian_lin_approx(double L)
{
	double Labs = fabs(L);

	if (Labs < 1.0)
	{
		return -0.375 * Labs + 0.6825;
	}
	else if ((Labs >= 1.0) && (Labs < 2.625))
	{
		return -0.1875 * Labs + 0.5;
	}
	else
	{
		return 0;
	}
}

int sign(double a)
{
	return (a <= 0) ? -1 : 1;
}
} // namespace ldpc