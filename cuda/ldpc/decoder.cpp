#include "../ldpcsim.h"

using namespace ldpc;

/*
*	Decoder device
*/
//Init constructor
__host__ ldpc_decoder_device::ldpc_decoder_device(cuda_ptr<ldpc_code_device> &pCode, size_t pI, bool pEarlyTerm)
    : mLdpcCode(pCode), mMaxIter(pI), mEarlyTerm(pEarlyTerm)
    , mLv2c(pCode->layers().size() * pCode->nnz()), mLc2v(pCode->layers().size() * pCode->nnz())
    , mLc2vPre(pCode->layers().size() * pCode->nnz()), mLSum(pCode->nnz())
    , mF(pCode->max_dc()), mB(pCode->max_dc())
    , mLLRIn(pCode->nc()), mLLROut(pCode->nc())
    , mSynd(pCode->mc()), mCO(pCode->nc())
    , mIter(0), mIsCW(false), FBREF(nullptr)
{
    cudaMallocManaged(&FBREF, mLdpcCode->max_dc());

    cudaDeviceSynchronize();
    int dev = -1;
    cudaGetDevice(&dev);
    cudaMemPrefetchAsync(FBREF, mLdpcCode->max_dc(), dev, NULL);

    //mem_prefetch();
}

//copy constructor
__host__ ldpc_decoder_device::ldpc_decoder_device(const ldpc_decoder_device &pCopy)
    : mLdpcCode(pCopy.mLdpcCode), mMaxIter(pCopy.mMaxIter), mEarlyTerm(pCopy.mEarlyTerm)
    , mLv2c(pCopy.mLv2c), mLc2v(pCopy.mLc2v), mLSum(pCopy.mLSum)
    , mLc2vPre(pCopy.mLc2vPre)
    , mF(pCopy.mF), mB(pCopy.mB)
    , mLLRIn(pCopy.mLLRIn), mLLROut(pCopy.mLLROut)
    , mSynd(pCopy.mSynd), mCO(pCopy.mCO)
    , mIter(0), mIsCW(false), FBREF(nullptr)
{
    cudaMallocManaged(&FBREF, mLdpcCode->max_dc());

    cudaDeviceSynchronize();
    int dev = -1;
    cudaGetDevice(&dev);
    cudaMemPrefetchAsync(FBREF, mLdpcCode->max_dc(), dev, NULL);

    //mem_prefetch();
}

//destructor
__host__ ldpc_decoder_device::~ldpc_decoder_device()
{
    if (FBREF != nullptr)
    {
        cudaFree(FBREF);
    }
}

//copy/move assignment operator
__host__ ldpc_decoder_device &ldpc_decoder_device::operator=(ldpc_decoder_device pCopy) noexcept
{
    swap(mLdpcCode, pCopy.mLdpcCode);
    swap(mLv2c, pCopy.mLv2c);
    swap(mLc2v, pCopy.mLc2v);
    swap(mLSum, pCopy.mLSum);
    swap(mLc2vPre, pCopy.mLc2vPre);
    swap(mLLRIn, pCopy.mLLRIn);
    swap(mLLROut, pCopy.mLLROut);
    swap(mF, pCopy.mF);
    swap(mB, pCopy.mB);
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

//legacy cpu decoder
__host__ __device__ size_t ldpc_decoder_device::decode_legacy()
{
    size_t it;

    size_t *vn;
    size_t *cn;

    size_t vw;
    size_t cw;

    //initialize with llrs
    for (size_t i = 0; i < mLdpcCode->nnz(); i++)
    {
        mLv2c[i] = mLLRIn[mLdpcCode->c()[i]];
    }

    it = 0;
    while (it < mMaxIter)
    {
        for (size_t i = 0; i < mLdpcCode->mc(); i++)
        {
            cw = mLdpcCode->cn()[i].size();
            cn = mLdpcCode->cn()[i].data();
            mF[0] = mLv2c[*cn];
            mB[cw - 1] = mLv2c[*(cn + cw - 1)];
            for (size_t j = 1; j < cw; j++)
            {
                mF[j] = jacobian(mF[j - 1], mLv2c[*(cn + j)]);
                mB[cw - 1 - j] = jacobian(mB[cw - j], mLv2c[*(cn + cw - j - 1)]);
            }

            mLc2v[*cn] = mB[1];
            mLc2v[*(cn + cw - 1)] = mF[cw - 2];
            for (size_t j = 1; j < cw - 1; j++)
            {
                mLc2v[*(cn + j)] = jacobian(mF[j - 1], mB[j + 1]);
            }
        }

        // VN node processing
        for (size_t i = 0; i < mLdpcCode->nc(); i++)
        {
            double tmp = mLLRIn[i];
            vw = mLdpcCode->vn()[i].size();
            vn = mLdpcCode->vn()[i].data();
            while (vw--)
            {
                tmp += mLc2v[*vn++];
            }
            vn = mLdpcCode->vn()[i].data();
            vw = mLdpcCode->vn()[i].size();
            while (vw--)
            {
                mLv2c[*vn] = tmp - mLc2v[*vn];
                vn++;
            }
        }

        // app calculation
        for (size_t i = 0; i < mLdpcCode->nc(); i++)
        {
            mLLROut[i] = mLLRIn[i];
            vn = mLdpcCode->vn()[i].data();
            vw = mLdpcCode->vn()[i].size();
            while (vw--)
            {
                mLLROut[i] += mLc2v[*vn++];
            }
            mCO[i] = (mLLROut[i] <= 0);
        }

        it++;

        if (mEarlyTerm)
        {
            if (is_codeword_legacy())
            {
                break;
            }
        }
    }

    return it;
}

__host__ __device__ bool ldpc_decoder_device::is_codeword_legacy()
{
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

    return true;
}
