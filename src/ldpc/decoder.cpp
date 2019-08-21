#include "../sim/ldpcsim.h"

namespace ldpc
{
/*
*	Decoder device
*/
//Init constructor
__host__ ldpc_decoder::ldpc_decoder(cuda_ptr<ldpc_code_device> &pCode, std::size_t pI, bool pEarlyTerm)
    : mLdpcCode(pCode), mMaxIter(pI), mEarlyTerm(pEarlyTerm),
      mLv2c(pCode->nnz()), mLc2v(pCode->nnz()),
      mExMsgCN(pCode->max_dc()),
      mLLRIn(pCode->nc()), mLLROut(pCode->nc()),
      mSynd(pCode->mc()), mCO(pCode->nc()),
      mIter(0), mIsCW(false)
{
}

//prefetch memory
__host__ void ldpc_decoder::mem_prefetch()
{
    mLc2v.mem_prefetch();
    mLv2c.mem_prefetch();

    mExMsgCN.mem_prefetch();

    mLLRIn.mem_prefetch();
    mLLROut.mem_prefetch();
    mCO.mem_prefetch();

    mSynd.mem_prefetch();
}

//calc syndrome & check if codeword
//calls global function
__host__ __device__ bool ldpc_decoder::is_codeword()
{
    mIsCW = true;

    //calc syndrome
    cudakernel::decoder::calc_synd<<<get_num_size(mLdpcCode->mc(), NUMK_THREADS), NUMK_THREADS>>>(this);
    cudaDeviceSynchronize();

    return mIsCW;
}

} // namespace ldpc