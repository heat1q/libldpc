#include "ldpc.h"

#include <math.h>

using namespace ldpc;


ldpc_code::ldpc_code(const char* pFileName, const char* pClFile, const bool pMgd)
	: cuda_mgd(pMgd)
{
	if (mIsMgd)
	{
		setup_mgd(pFileName);
		setup_layers_mgd(pClFile);
		prefetch();
	}
	else
	{
		setup(pFileName);
		setup_layers(pClFile);
	}
}


ldpc_code::~ldpc_code()
{
	if (mIsMgd)
	{
		destroy_mgd();
	}
	else
	{
		destroy();
	}
}


void ldpc_code::setup(const char* pFileName)
{
	//init
	mPuncture = nullptr;
	mShorten = nullptr;
	mCW = nullptr;
	mCN = nullptr;
	mVW = nullptr;
	mVN = nullptr;
	mR = nullptr;
	mC = nullptr;
	mLW = nullptr;
	mLayers = nullptr;

	try
	{
		FILE *fp;

		fp = fopen(pFileName, "r");
		if(!fp)
			throw std::runtime_error("can not open codefile for reading.");

		fscanf(fp, "nc: %lu\n", &mN);
		fscanf(fp, "mc: %lu\n", &mM);
		fscanf(fp, "nct: %lu\n", &mNCT);
		fscanf(fp, "mct: %lu\n", &mMCT);
		fscanf(fp,  "nnz: %lu\n", &mNNZ);
		mK = mN-mM;
		mKCT = mNCT-mMCT;

		fscanf(fp, "puncture [%lu]: ", &(mNumPuncture));
		mNumPunctureSys = 0;
		mNumPuncturePar = 0;
		if(mNumPuncture != 0)
		{
			mPuncture = new size_t[mNumPuncture];
			for(size_t i = 0; i < mNumPuncture; i++)
			{
				fscanf(fp, " %lu ", &(mPuncture[i]));
				if(mPuncture[i] < mK)
					mNumPunctureSys++;
				else
					mNumPuncturePar++;
			}
		}

		fscanf(fp, "shorten [%lu]: ", &mNumShorten);
		if(mNumShorten != 0)
		{
			mShorten = new size_t[mNumShorten];
			for(size_t i = 0; i < mNumShorten; i++)
				fscanf(fp, " %lu ", &(mShorten[i]));
		}


		size_t* cwTmp;
		size_t* vwTmp;
		mCW = new uint64_t[mM] ();
		cwTmp = new uint64_t[mM] ();
		mVW = new uint64_t[mN] ();
		vwTmp = new uint64_t[mN] ();
		mR = new uint64_t[mNNZ] ();
		mC = new uint64_t[mNNZ] ();


		for(size_t i = 0; i < mNNZ; i++)
		{
			fscanf(fp, "%lu %lu\n", &(mR[i]), &(mC[i]));
			mCW[mR[i]]++;
			mVW[mC[i]]++;
		}

		mCN = new size_t*[mM] ();
		for(size_t i = 0; i < mM; i++)
			mCN[i] = new size_t[mCW[i]] ();

		mVN = new size_t*[mN] ();
		for(size_t i = 0; i < mN; i++)
			mVN[i] = new size_t[mVW[i]] ();

		for(size_t i = 0; i < mNNZ; i++)
		{
			mCN[mR[i]][cwTmp[mR[i]]++] = i;
			mVN[mC[i]][vwTmp[mC[i]]++] = i;
		}

		delete[] cwTmp;
		delete[] vwTmp;

		mMaxDC = 0;
		for(size_t i = 0; i < mM; i++)
		{
			if(mCW[i] > mMaxDC)
				mMaxDC = mCW[i];
		}

		fclose(fp);
	}
	catch(std::exception &e)
	{
		std::cout << "Error: " << e.what() << "\n";
		destroy();

		exit(EXIT_FAILURE);
	}
}


void ldpc_code::setup_layers(const char* pClFile)
{
    FILE *fp = fopen(pClFile, "r");
    if(!fp)
        throw std::runtime_error("Can not open layer file");

    fscanf(fp, "nl: %lu\n", &mNL);

    mLW = new uint64_t[mNL];
    mLayers = new uint64_t*[mNL];

    for (size_t i = 0; i < mNL; ++i)
    {
        fscanf(fp, "cn[i]: %lu\n", &(mLW[i]));
        mLayers[i] = new uint64_t[mLW[i]];
        for (size_t j = 0; j < mLW[i]; ++j)
            fscanf(fp, "%lu\n", &(mLayers[i][j]));
    }
    fclose(fp);
}


void ldpc_code::setup_mgd(const char* pFileName)
{
	//init
    mPuncture = nullptr;
    mShorten = nullptr;
    mCW = nullptr;
    mCN = nullptr;
    mVW = nullptr;
    mVN = nullptr;
    mR = nullptr;
    mC = nullptr;
    mLW = nullptr;
    mLayers = nullptr;

    try
    {
        FILE *fp;

        fp = fopen(pFileName, "r");
        if(!fp)
            throw std::runtime_error("can not open codefile for reading.");

        fscanf(fp, "nc: %lu\n", &mN);
        fscanf(fp, "mc: %lu\n", &mM);
        fscanf(fp, "nct: %lu\n", &mNCT);
        fscanf(fp, "mct: %lu\n", &mMCT);
        fscanf(fp,  "nnz: %lu\n", &mNNZ);
        mK = mN-mM;
        mKCT = mNCT-mMCT;

        fscanf(fp, "puncture [%lu]: ", &(mNumPuncture));
        mNumPunctureSys = 0;
        mNumPuncturePar = 0;
        if(mNumPuncture != 0)
        {
			cudaMallocManaged(&mPuncture, mNumPuncture*sizeof(size_t));
            for(size_t i = 0; i < mNumPuncture; i++)
            {
                fscanf(fp, " %lu ", &(mPuncture[i]));
                if(mPuncture[i] < mK)
                    mNumPunctureSys++;
                else
                    mNumPuncturePar++;
            }
        }

        fscanf(fp, "shorten [%lu]: ", &mNumShorten);
        if(mNumShorten != 0)
        {
			cudaMallocManaged(&mShorten, mNumShorten*sizeof(size_t));
            for(size_t i = 0; i < mNumShorten; i++)
                fscanf(fp, " %lu ", &(mShorten[i]));
        }

        size_t* cwTmp;
        size_t* vwTmp;

		cwTmp = new size_t[mM]();
		vwTmp = new size_t[mN]();

		cudaMallocManaged(&mCW, mM*sizeof(size_t));
		cudaMallocManaged(&mVW, mN*sizeof(size_t));
		cudaMallocManaged(&mR, mNNZ*sizeof(size_t));
		cudaMallocManaged(&mC, mNNZ*sizeof(size_t));

		for(size_t i = 0; i < mM; i++)
		{
	        mCW[i] = 0;
	        cwTmp[i] = 0;
	    }
	    for(size_t i = 0; i < mN; i++)
		{
	        mVW[i] = 0;
	        vwTmp[i] = 0;
	    }

        for(size_t i = 0; i < mNNZ; i++)
        {
            fscanf(fp, "%lu %lu\n", &(mR[i]), &(mC[i]));
            mCW[mR[i]]++;
            mVW[mC[i]]++;
        }

		cudaMallocManaged(&mCN, mM*sizeof(size_t*));
        for(size_t i = 0; i < mM; i++)
			cudaMallocManaged(&mCN[i], mCW[i]*sizeof(size_t));

		cudaMallocManaged(&mVN, mN*sizeof(size_t*));
        for(size_t i = 0; i < mN; i++)
			cudaMallocManaged(&mVN[i], mVW[i]*sizeof(size_t));

        for(size_t i = 0; i < mNNZ; i++)
        {
            mCN[mR[i]][cwTmp[mR[i]]++] = i;
            mVN[mC[i]][vwTmp[mC[i]]++] = i;
        }

		delete[] cwTmp;
		delete[] vwTmp;

        mMaxDC = 0;
        for(size_t i = 0; i < mM; i++)
        {
            if(mCW[i] > mMaxDC)
                mMaxDC = mCW[i];
        }

        fclose(fp);
    }
    catch(std::exception &e)
    {
        std::cout << "Error: " << e.what() << "\n";
        destroy_mgd();

        exit(EXIT_FAILURE);
    }
}


void ldpc_code::setup_layers_mgd(const char* pClFile)
{
    FILE *fp = fopen(pClFile, "r");
    if(!fp)
        throw std::runtime_error("Can not open layer file");

    fscanf(fp, "nl: %lu\n", &mNL);

	cudaMallocManaged(&mLW, mNL*sizeof(uint64_t));
	cudaMallocManaged(&mLayers, mNL*sizeof(uint64_t*));

    for (size_t i = 0; i < mNL; ++i)
    {
        fscanf(fp, "cn[i]: %lu\n", &(mLW[i]));
		cudaMallocManaged(&mLayers[i], mLW[i]*sizeof(uint64_t));
        for (size_t j = 0; j < mLW[i]; ++j)
            fscanf(fp, "%lu\n", &(mLayers[i][j]));
    }
    fclose(fp);
}


void ldpc_code::prefetch()
{
	cudaDeviceSynchronize();

	int dev = -1;
	cudaGetDevice(&dev);

	if(mNumPuncture != 0)	{ cudaMemPrefetchAsync(mPuncture, sizeof(size_t)*mNumPuncture, dev, NULL); }
	if(mNumShorten != 0) { cudaMemPrefetchAsync(mShorten, sizeof(size_t)*mNumShorten, dev, NULL); }


	for(size_t i = 0; i < mM; i++) {
		cudaMemPrefetchAsync(mCN[i], sizeof(size_t)*mCW[i], dev, NULL);
	}
	cudaMemPrefetchAsync(mCN, sizeof(size_t*)*mM, dev, NULL);


	for(size_t i = 0; i < mM; i++) {
		cudaMemPrefetchAsync(mVN[i], sizeof(size_t)*mVW[i], dev, NULL);
	}
	cudaMemPrefetchAsync(mVN, sizeof(size_t*)*mN, dev, NULL);


	for (size_t i = 0; i < mNL; ++i) {
		cudaMemPrefetchAsync(mLayers[i], sizeof(uint64_t)*mLW[i], dev, NULL);
	}
	cudaMemPrefetchAsync(mLayers, sizeof(uint64_t*)*mNL, dev, NULL);
	cudaMemPrefetchAsync(mLW, sizeof(uint64_t)*mNL, dev, NULL);

	cudaMemPrefetchAsync(mCW, sizeof(size_t)*mM, dev, NULL);
	cudaMemPrefetchAsync(mVW, sizeof(size_t)*mN, dev, NULL);
	cudaMemPrefetchAsync(mR, sizeof(size_t)*mNNZ, dev, NULL);
	cudaMemPrefetchAsync(mC, sizeof(size_t)*mNNZ, dev, NULL);

	cudaMemPrefetchAsync(this, sizeof(ldpc_code), dev, NULL);
}


void ldpc_code::destroy()
{
    if (mVN != nullptr)
    {
        for(size_t i = 0; i < mN; i++) { delete[] mVN[i]; }
        delete[] mVN;
    }

    if (mVN != nullptr)
    {
        for(size_t i = 0; i < mM; i++) { delete[] mCN[i]; }
        delete[] mCN;
    }

    if (mVW != nullptr) { delete[] mVW; }
    if (mCW != nullptr) { delete[] mCW; }
    if (mR != nullptr) { delete[] mR; }
    if (mC != nullptr) { delete[] mC; }
    if (mPuncture != nullptr) { delete[] mPuncture; }
    if (mShorten != nullptr) { delete[] mShorten; }
    if (mLayers != nullptr)
    {
        for(size_t i = 0; i < mNL; i++) { delete[] mLayers[i]; }
        delete[] mLayers;
    }
    if (mLW != nullptr) { delete[] mLW; }
}


void ldpc_code::destroy_mgd()
{
    if (mVN != nullptr)
    {
        for(size_t i = 0; i < mN; i++) { cudaFree(mVN[i]); }
        cudaFree(mVN);
    }
    if (mVN != nullptr)
    {
        for(size_t i = 0; i < mM; i++) { cudaFree(mCN[i]); }
        cudaFree(mCN);
    }
    if (mVW != nullptr) { cudaFree(mVW); }
    if (mCW != nullptr) { cudaFree(mCW); }
    if (mR != nullptr) { cudaFree(mR); }
    if (mC != nullptr) { cudaFree(mC); }
    if (mPuncture != nullptr) { cudaFree(mPuncture); }
    if (mShorten != nullptr) { cudaFree(mShorten); }
    if (mLayers != nullptr)
	{
        for(size_t i = 0; i < mNL; i++) { cudaFree(mLayers[i]); }
        cudaFree(mLayers);
    }
    if (mLW != nullptr) { cudaFree(mLW); }
}

void ldpc_code::print()
{
    std::cout << "=========== LDPC ===========" << "\n";
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
    std::cout << "=========== LDPC: END ===========" << "\n";

    printf("=========== LDPC LAYERS ===========\n");
    printf("nl: %lu\n", mNL);
    for (size_t i = 0; i < mNL; ++i)
    {
        printf("cn[%lu]: %lu\n", i, mLW[i]);
        //printVector<uint64_t>(mLayers[i], mLW[i]);
        printf("\n");
    }
    printf("========= LDPC LAYERS: END ========\n");
}

__host__ __device__ uint64_t ldpc_code::nc() const { return mN; }
__host__ __device__ uint64_t ldpc_code::kc() const { return mK; }
__host__ __device__ uint64_t ldpc_code::mc() const { return mM; }
__host__ __device__ uint64_t ldpc_code::nnz() const { return mNNZ; }
__host__ __device__ size_t *ldpc_code::cw() const { return mCW; }
__host__ __device__ size_t *ldpc_code::vw() const { return mVW; }
__host__ __device__ size_t **ldpc_code::cn() const { return mCN; }
__host__ __device__ size_t **ldpc_code::vn() const { return mVN; }
__host__ __device__ size_t *ldpc_code::r() const { return mR; }
__host__ __device__ size_t *ldpc_code::c() const { return mC; }
__host__ __device__ uint64_t ldpc_code::nct() const { return mNCT; }
__host__ __device__ uint64_t ldpc_code::mct() const { return mMCT; }
__host__ __device__ size_t *ldpc_code::puncture() const { return mPuncture; }
__host__ __device__ size_t ldpc_code::num_puncture() const { return mNumPuncture; }
__host__ __device__ size_t *ldpc_code::shorten() const { return mShorten; }
__host__ __device__ size_t ldpc_code::num_shorten() const { return mNumShorten; }
__host__ __device__ uint64_t ldpc_code::kct() const { return mKCT; }
__host__ __device__ size_t ldpc_code::max_dc() const { return mMaxDC; }
__host__ __device__ uint64_t ldpc_code::nl() const { return mNL; }
__host__ __device__ uint64_t *ldpc_code::lw() const { return mLW; }
__host__ __device__ uint64_t **ldpc_code::layers() const { return mLayers; }


void ldpc::dec2bin(uint64_t val, uint8_t m)
{
    for(size_t i = 0; i < m; i++)
        printf("%lu", (val>>(m-i-1) & 0x01));
}

__host__ __device__ double ldpc::jacobian(const double& L1, const double& L2)
{
#ifdef CN_APPROX_LIN
    return sign(L1) * sign(L2) * fmin(fabs(L1),fabs(L2)) + jacobian_lin_approx(L1+L2) - jacobian_lin_approx(L1-L2);
#elif CN_APPROX_MINSUM
    return sign(L1) * sign(L2) * fmin(fabs(L1), fabs(L2));
#else
    return sign(L1) * sign(L2) * fmin(fabs(L1),fabs(L2)) + log((1+exp(-fabs(L1+L2)))/(1+exp(-fabs(L1-L2))));
#endif
}

__host__ __device__ double ldpc::jacobian_lin_approx(const double& L)
{
    double Labs = fabs(L);

    if(Labs < 1.0) {
        return -0.375 * Labs  + 0.6825;
    } else if((Labs >= 1.0) && (Labs < 2.625)) {
        return -0.1875 * Labs + 0.5;
    } else {
        return 0;
    }
}

__host__ __device__ int8_t ldpc::sign(const double& a)
{
    return (a <= 0) ? -1 : 1;
}


/*
* Code device
*/
//init constructor
ldpc_code_device::ldpc_code_device(const char* pFileName, const char* pClFile)
: mMaxDC(0)
{
	try
	{
		FILE *fpCode = fopen(pFileName, "r");
		if(!fpCode) { throw std::runtime_error("can not open codefile for reading."); }
		FILE *fpLayer = fopen(pClFile, "r");
		if(!fpLayer) { throw std::runtime_error("Can not open layer file"); }

		fscanf(fpCode, "nc: %lu\n", &mN);
		fscanf(fpCode, "mc: %lu\n", &mM);
		fscanf(fpCode, "nct: %lu\n", &mNCT);
		fscanf(fpCode, "mct: %lu\n", &mMCT);
		fscanf(fpCode,  "nnz: %lu\n", &mNNZ);
		mK = mN-mM;
		mKCT = mNCT-mMCT;

		fscanf(fpCode, "puncture [%lu]: ", &(mNumPuncture));
		mNumPunctureSys = 0;
		mNumPuncturePar = 0;
		if(mNumPuncture != 0)
		{
			mPuncture = vec_size_t(mNumPuncture);
			for(size_t i = 0; i < mNumPuncture; i++)
			{
				fscanf(fpCode, " %lu ", &(mPuncture[i]));
				if (mPuncture[i] < mK) { mNumPunctureSys++; }
				else { mNumPuncturePar++; }
			}
		}

		fscanf(fpCode, "shorten [%lu]: ", &mNumShorten);
		if(mNumShorten != 0)
		{
			mShorten = vec_size_t(mNumShorten);
			for(size_t i = 0; i < mNumShorten; i++)
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


		for(size_t i = 0; i < mNNZ; i++)
		{
			fscanf(fpCode, "%lu %lu\n", &(mR[i]), &(mC[i]));
			mCW[mR[i]]++;
			mVW[mC[i]]++;
		}

		mCN = mat_size_t(mM, vec_size_t());
		for(size_t i = 0; i < mM; i++)
		{
			mCN[i] = vec_size_t(mCW[i]);
		}

		mVN = mat_size_t(mN, vec_size_t());
		for(size_t i = 0; i < mN; i++)
		{
			mVN[i] = vec_size_t(mVW[i]);
		}

		for(size_t i = 0; i < mNNZ; i++)
		{
			mCN[mR[i]][cwTmp[mR[i]]++] = i;
			mVN[mC[i]][vwTmp[mC[i]]++] = i;
		}

		for(size_t i = 0; i < mM; i++)
		{
			if(mCW[i] > mMaxDC) { mMaxDC = mCW[i]; }
		}

		//setup layers
		fscanf(fpLayer, "nl: %lu\n", &mNL);

		mLW = vec_size_t(mNL);
		mLayers = mat_size_t(mNL, vec_size_t());

		for (size_t i = 0; i < mNL; ++i)
		{
			fscanf(fpLayer, "cn[i]: %lu\n", &(mLW[i]));
			mLayers[i] = vec_size_t(mLW[i]);
			for (size_t j = 0; j < mLW[i]; ++j)
				fscanf(fpLayer, "%lu\n", &(mLayers[i][j]));
		}

		fclose(fpLayer);
		fclose(fpCode);
	}
	catch(std::exception &e)
	{
		std::cout << "Error: " << e.what() << "\n";
		exit(EXIT_FAILURE);
	}
}


void ldpc_code_device::print()
{
    std::cout << "=========== LDPC ===========" << "\n";
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
    std::cout << "=========== LDPC: END ===========" << std::endl;

    printf("=========== LDPC LAYERS ===========\n");
    printf("nl: %lu\n", mNL);
    for (size_t i = 0; i < mNL; ++i)
    {
        printf("cn[%lu]: %lu\n", i, mLW[i]);
		for (auto x : mLayers[i]) {
			printf("%lu ", x);
		}
        printf("\n");
    }
    printf("========= LDPC LAYERS: END ========\n");
}


void ldpc_code_device::prefetch()
{
	//TODO
}
