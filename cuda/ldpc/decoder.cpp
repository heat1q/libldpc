#include "ldpc.h"

using namespace ldpc;


ldpc_decoder::ldpc_decoder(ldpc_code* pCode, const uint16_t pI, const bool pEarlyTerm)
	: mLdpcCode(pCode), mMaxIter(pI), mEarlyTerm(pEarlyTerm), cuda_mgd(true)
{
	setup();
	prefetch();
}

ldpc_decoder::~ldpc_decoder() { destroy(); }


void ldpc_decoder::setup()
{
    mLc2v = nullptr;
    mLv2c = nullptr;
    mF = nullptr;
    mB = nullptr;
    FBREF = nullptr;
    mLSum = nullptr;
    mLc2vPre = nullptr;
    mCO = nullptr;
    mSynd = nullptr;
    mLLRIn = nullptr;
    mLLROut = nullptr;

    const uint64_t numLayers = mLdpcCode->nl();

    try
    {
        //num layers times num nnz
        cudaMallocManaged(&mLc2v, sizeof(double)*numLayers*mLdpcCode->nnz());
        if (mLc2v == NULL || mLc2v == nullptr) { throw std::runtime_error("mLc2v alloc failed."); }
        cudaMallocManaged(&mLv2c, sizeof(double)*numLayers*mLdpcCode->nnz());
        if (mLv2c == NULL || mLv2c == nullptr) { throw std::runtime_error("mLv2c alloc failed."); }
        cudaMallocManaged(&mLc2vPre, sizeof(double)*numLayers*mLdpcCode->nnz());
        if (mLc2vPre == NULL || mLc2vPre == nullptr) { throw std::runtime_error("mLc2vPre alloc failed."); }

        cudaMallocManaged(&mF, sizeof(double)*numLayers*mLdpcCode->max_dc());
        if (mF == NULL || mF == nullptr) { throw std::runtime_error("mF alloc failed."); }
        cudaMallocManaged(&mB, sizeof(double)*numLayers*mLdpcCode->max_dc());
        if (mB == NULL || mB == nullptr) { throw std::runtime_error("mB alloc failed."); }

        cudaMallocManaged(&FBREF, mLdpcCode->max_dc());
        if (FBREF == NULL || FBREF == nullptr) { throw std::runtime_error("FBREF alloc failed."); }

        cudaMallocManaged(&mLSum, sizeof(double)*mLdpcCode->nnz());
        if (mLSum == NULL || mLSum == nullptr) { throw std::runtime_error("mLSum alloc failed."); }

        cudaMallocManaged(&mLLRIn, sizeof(double)*mLdpcCode->nc());
        if (mLLRIn == NULL || mLLRIn == nullptr) { throw std::runtime_error("mLLRIn alloc failed."); }
        cudaMallocManaged(&mLLROut, sizeof(double)*mLdpcCode->nc());
        if (mLLROut == NULL || mLLROut == nullptr) { throw std::runtime_error("mLLROut alloc failed."); }
        cudaMallocManaged(&mCO, sizeof(bits_t)*mLdpcCode->nc());
        if (mCO == NULL || mCO == nullptr) { throw std::runtime_error("mCO alloc failed."); }

        cudaMallocManaged(&mSynd, sizeof(bits_t)*mLdpcCode->mc());
        if (mSynd == NULL || mSynd == nullptr) { throw std::runtime_error("mSynd alloc failed."); }
    }
    catch (std::exception& e)
    {
        std::cout << "Error: " << e.what() << "\n";
        destroy();
        exit(EXIT_FAILURE);
    }
}

void ldpc_decoder::prefetch()
{
    cudaDeviceSynchronize();

    int dev = -1;
    cudaGetDevice(&dev);

    const uint64_t numLayers = mLdpcCode->nl();

    cudaMemPrefetchAsync(mLc2v, sizeof(double)*numLayers*mLdpcCode->nnz(), dev, NULL);
    cudaMemPrefetchAsync(mLv2c, sizeof(double)*numLayers*mLdpcCode->nnz(), dev, NULL);
    cudaMemPrefetchAsync(mLc2vPre, sizeof(double)*numLayers*mLdpcCode->nnz(), dev, NULL);

    cudaMemPrefetchAsync(FBREF, mLdpcCode->max_dc(), dev, NULL);

    cudaMemPrefetchAsync(mLSum, sizeof(double)*mLdpcCode->nnz(), dev, NULL);

    cudaMemPrefetchAsync(mLLRIn, sizeof(double)*mLdpcCode->nc(), dev, NULL);
    cudaMemPrefetchAsync(mLLROut, sizeof(double)*mLdpcCode->nc(), dev, NULL);
    cudaMemPrefetchAsync(mCO, sizeof(double)*mLdpcCode->nc(), dev, NULL);

    cudaMemPrefetchAsync(mSynd, sizeof(double)*mLdpcCode->mc(), dev, NULL);

    cudaMemPrefetchAsync(this, sizeof(ldpc_decoder), dev, NULL);
}

void ldpc_decoder::destroy()
{
    if (mLc2v != nullptr) { cudaFree(mLc2v); }
    if (mLv2c != nullptr) { cudaFree(mLv2c); }
    if (mLc2vPre != nullptr) { cudaFree(mLc2vPre); }
    if (mF != nullptr) { cudaFree(mF); }
    if (mB != nullptr) { cudaFree(mB); }
    if (FBREF != nullptr) { cudaFree(FBREF); }
    if (mLSum != nullptr) { cudaFree(mLSum); }
    if (mCO != nullptr) { cudaFree(mCO); }
    if (mSynd != nullptr) { cudaFree(mSynd); }
    if (mLLRIn != nullptr) { cudaFree(mLLRIn); }
    if (mLLROut != nullptr) { cudaFree(mLLROut); }
}


__host__ __device__ bool ldpc_decoder::is_codeword()
{
    mIsCW = true;

    //calc syndrome
    //cudakernel::decoder::calc_synd<<<get_num_size(mLdpcCode->mc(), NUM_THREADS), NUM_THREADS>>>(this);
    cudaDeviceSynchronize();

    return mIsCW;
}


__host__ __device__ bool ldpc_decoder::is_codeword_legacy()
{
    bool is_codeword = true;

    //calc syndrome
    bits_t s;
    for (size_t i = 0; i < mLdpcCode->mc(); i++)
    {
        s = 0;
        for (size_t j = 0; j < mLdpcCode->cw()[i]; j++)
            s ^= mCO[mLdpcCode->c()[mLdpcCode->cn()[i][j]]];

        if (s)
        {
            return false;
        }
    }

    return is_codeword;
}


uint16_t ldpc_decoder::decode_legacy()
{
    uint16_t it;

    size_t* vn;
    size_t* cn;

    size_t vw;
    size_t cw;

    //initialize with llrs
    for(size_t i = 0; i < mLdpcCode->nnz(); i++) {
        mLv2c[i] = mLLRIn[mLdpcCode->c()[i]];
    }

    it = 0;
    while(it < mMaxIter) {
        for(size_t i = 0; i < mLdpcCode->mc(); i++) {
            cw = mLdpcCode->cw()[i];
            cn = mLdpcCode->cn()[i];
            mF[0] = mLv2c[*cn];
            mB[cw-1] = mLv2c[*(cn+cw-1)];
            for(size_t j = 1; j < cw; j++) {
                mF[j] = jacobian(mF[j-1], mLv2c[*(cn+j)]);
                mB[cw-1-j] = jacobian(mB[cw-j], mLv2c[*(cn + cw-j-1)]);
            }

            mLc2v[*cn] = mB[1];
            mLc2v[*(cn+cw-1)] = mF[cw-2];
            for(size_t j = 1; j < cw-1; j++) {
                mLc2v[*(cn+j)] = jacobian(mF[j-1], mB[j+1]);
            }
        }

        // VN node processing
        for(size_t i = 0; i < mLdpcCode->nc(); i++) {
            double tmp = mLLRIn[i];
            vw = mLdpcCode->vw()[i];
            vn = mLdpcCode->vn()[i];
            while(vw--) {
                tmp += mLc2v[*vn++];
            }
            vn = mLdpcCode->vn()[i];
            vw = mLdpcCode->vw()[i];
            while(vw--) {
                mLv2c[*vn] = tmp - mLc2v[*vn];
                vn++;
            }
        }

        // app calculation
        for(size_t i = 0; i < mLdpcCode->nc(); i++) {
            mLLROut[i] = mLLRIn[i];
            vn = mLdpcCode->vn()[i];
            vw = mLdpcCode->vw()[i];
            while(vw--) {
                mLLROut[i] += mLc2v[*vn++];
            }
            mCO[i] = (mLLROut[i] <= 0);
        }

        it++;

        if (mEarlyTerm) {
            if (is_codeword_legacy()) {
                break;
            }
        }
    }

    return it;
}


uint16_t ldpc_decoder::decode_layered_legacy()
{
    size_t* vn;
    size_t* cn;

    size_t vw;
    size_t cw;

    size_t iNNZ;
    size_t iDC;

    //initialize
    for (size_t i = 0; i < mLdpcCode->nnz(); ++i)
    {
        mLSum[i] = 0.0;
        for (size_t l = 0; l < mLdpcCode->nl(); ++l)
        {
            mLc2v[l*mLdpcCode->nnz()+i] = 0.0;
            mLv2c[l*mLdpcCode->nnz()+i] = 0.0;
            mLc2vPre[l*mLdpcCode->nnz()+i] = 0.0;
        }
    }

    uint16_t I = 0;
    while (I < mMaxIter)
    {
        for (size_t l = 0; l < mLdpcCode->nl(); ++l)
        {
            iNNZ = l*mLdpcCode->nnz();
            iDC = l*mLdpcCode->max_dc();

            // VN node intialization
            for(size_t i = 0; i < mLdpcCode->nc(); i++)
            {
                double tmp = mLLRIn[i];
                vw = mLdpcCode->vw()[i];
                vn = mLdpcCode->vn()[i];
                while(vw--)
                    tmp += mLSum[*vn++];

                vn = mLdpcCode->vn()[i];
                vw = mLdpcCode->vw()[i];
                while(vw--)
                {
                    mLv2c[iNNZ + *vn] = tmp - mLc2v[iNNZ + *vn];
                    ++vn;
                }
            }

            //CN processing
            for(size_t i = 0; i < mLdpcCode->lw()[l]; i++)
            {
                cw = mLdpcCode->cw()[mLdpcCode->layers()[l][i]];
                cn = mLdpcCode->cn()[mLdpcCode->layers()[l][i]];
                mF[iDC] = mLv2c[iNNZ + *cn];
                mB[iDC + cw-1] = mLv2c[iNNZ + *(cn+cw-1)];
                for(size_t j = 1; j < cw; j++)
                {
                    mF[iDC + j] = jacobian(mF[iDC + j-1], mLv2c[iNNZ + *(cn+j)]);
                    mB[iDC + cw-1-j] = jacobian(mB[iDC + cw-j], mLv2c[iNNZ + *(cn + cw-j-1)]);
                }

                mLc2v[iNNZ + *cn] = mB[iDC + 1];
                mLc2v[iNNZ + *(cn+cw-1)] = mF[iDC + cw-2];

                for(size_t j = 1; j < cw-1; j++)
                    mLc2v[iNNZ + *(cn+j)] = jacobian(mF[iDC + j-1], mB[iDC + j+1]);
            }

            //update the llr sum of layers, by replacing old llr of lyr l with new value
            for (size_t i = 0; i < mLdpcCode->nnz(); ++i)
            {
                mLSum[i] += mLc2v[iNNZ + i] - mLc2vPre[iNNZ + i];
                mLc2vPre[iNNZ + i] = mLc2v[iNNZ + i];
            }

            // app calculation
            for(size_t i = 0; i < mLdpcCode->nc(); ++i)
            {
                mLLROut[i] = mLLRIn[i];
                vn = mLdpcCode->vn()[i];
                vw = mLdpcCode->vw()[i];
                while(vw--)
                    mLLROut[i] += mLSum[*vn++];
                mCO[i] = (mLLROut[i] <= 0);
            }

            if (mEarlyTerm)
            {
                if (is_codeword_legacy())
                {
                    return I;
                }
            }
        }

        ++I;
    }

    return I;
}


uint16_t ldpc_decoder::decode_layered()
{
    //cudakernel::decoder::decode_layered<<<1, 1>>>(this);
    cudaDeviceSynchronize();

    return mIter;
}



/*
*	Decoder device
*/
//Init constructor
__host__ ldpc_decoder_device::ldpc_decoder_device(cudamgd_ptr<ldpc_code_device> pCode, size_t pI, bool pEarlyTerm)
: mLdpcCode(pCode), mMaxIter(pI), mEarlyTerm(pEarlyTerm)
, mLv2c(pCode->layers().size()*pCode->nnz())
, mLc2v(pCode->layers().size()*pCode->nnz())
, mLc2vPre(pCode->layers().size()*pCode->nnz())
, mLSum(pCode->nnz())
, mLLRIn(pCode->nc()), mLLROut(pCode->nc())
, mSynd(pCode->mc()), mCO(pCode->nc())
, mIter(0), mIsCW(false), FBREF(nullptr)
{
	cudaMallocManaged(&FBREF, mLdpcCode->max_dc());

	cudaDeviceSynchronize();
	int dev = -1;
	cudaGetDevice(&dev);
	cudaMemPrefetchAsync(FBREF, mLdpcCode->max_dc(), dev, NULL);

	mem_prefetch();
}

//copy constructor
__host__ ldpc_decoder_device::ldpc_decoder_device(const ldpc_decoder_device& pCopy)
: mLdpcCode(pCopy.mLdpcCode), mMaxIter(pCopy.mMaxIter), mEarlyTerm(pCopy.mEarlyTerm)
, mLv2c(pCopy.mLv2c), mLc2v(pCopy.mLc2v), mLSum(pCopy.mLSum), mLc2vPre(pCopy.mLc2vPre)
, mLLRIn(pCopy.mLLRIn), mLLROut(pCopy.mLLROut)
, mSynd(pCopy.mSynd), mCO(pCopy.mCO)
, mIter(0), mIsCW(false), FBREF(nullptr)
{
	cudaMallocManaged(&FBREF, mLdpcCode->max_dc());

	cudaDeviceSynchronize();
	int dev = -1;
	cudaGetDevice(&dev);
	cudaMemPrefetchAsync(FBREF, mLdpcCode->max_dc(), dev, NULL);

	mem_prefetch();
}

//destructor
__host__ ldpc_decoder_device::~ldpc_decoder_device()
{
	if (FBREF != nullptr) { cudaFree(FBREF); }
}

//copy/move assignment operator
__host__ ldpc_decoder_device& ldpc_decoder_device::operator=(ldpc_decoder_device pCopy) noexcept
{
	swap(mLdpcCode, pCopy.mLdpcCode);
	swap(mLv2c, pCopy.mLv2c);
	swap(mLc2v, pCopy.mLc2v);
	swap(mLSum, pCopy.mLSum);
	swap(mLc2vPre, pCopy.mLc2vPre);
	swap(mLLRIn, pCopy.mLLRIn);
	swap(mLLROut, pCopy.mLLROut);
	swap(mSynd, pCopy.mSynd);
	swap(mCO, pCopy.mCO);
	swap(mIter, pCopy.mIter);
	swap(FBREF, pCopy.FBREF);
	swap(mIsCW, pCopy.mIsCW);
	swap(mMaxIter, pCopy.mMaxIter);
	swap(mEarlyTerm, pCopy.mEarlyTerm);

	return *this;
}

//prefetch memory
__host__ void ldpc_decoder_device::mem_prefetch()
{
	mLc2v.mem_prefetch();
	mLv2c.mem_prefetch();
	mLc2vPre.mem_prefetch();

	mLSum.mem_prefetch();

	mLLRIn.mem_prefetch();
	mLLROut.mem_prefetch();
	mCO.mem_prefetch();

	mSynd.mem_prefetch();
}


//calc syndrome & check if codeword
//calls global function
__host__ __device__ bool ldpc_decoder_device::is_codeword()
{
    mIsCW = true;

    //calc syndrome
    cudakernel::decoder::calc_synd<<<get_num_size(mLdpcCode->mc(), NUM_THREADS), NUM_THREADS>>>(this);
    cudaDeviceSynchronize();

    return mIsCW;
}

//start decoding
//calls global function
__host__ __device__ size_t ldpc_decoder_device::decode_layered()
{
	cudakernel::decoder::decode_layered<<<1, 1>>>(this);
    cudaDeviceSynchronize();

    return mIter;
}
