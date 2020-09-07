#include "decoder.h"
#include "../sim/ldpcsim.h"

namespace ldpc
{
/**
 * @brief Construct a new ldpc decoder::ldpc decoder object
 * 
 * @param pCode 
 * @param pI 
 * @param pEarlyTerm 
 */
ldpc_decoder::ldpc_decoder(ldpc_code *pCode, ldpc_sim *pSim, unsigned pI, bool pEarlyTerm)
    : mLdpcCode(pCode), mSim(pSim), 
      mLv2c(pCode->nnz()), mLc2v(pCode->nnz()),
      mExMsgCN(pCode->max_dc()),
      mLLRIn(pCode->nc()), mLLROut(pCode->nc()),
      mSynd(pCode->mc()), mCO(pCode->nc()),
      mMaxIter(pI), mEarlyTerm(pEarlyTerm)
{
}

void ldpc_decoder::calc_llrs(const vec_double_t &y, double sigma2)
{
    //puncturing & shortening
    if (mLdpcCode->puncture().size() != 0)
    {
        for (auto p: mLdpcCode->puncture())
        {
            mLLRIn[p] = 1e-9 * (1 - 2*(rand() % 2)); // set random signed low values, i.e. simulate erasure LLR=0
        }
    }
    if (mLdpcCode->shorten().size() != 0)
    {
        for (auto s: mLdpcCode->shorten())
        {
            mLLRIn[s] = 1e9; // simulate certain bit
        }
    }

    //bpsk
    for (u64 i = 0; i < mLdpcCode->nct(); ++i)
    {
        mLLRIn[mSim->bits_pos()[i]] = 2 * y[i] / sigma2;
    }
}

/**
 * @brief Standard BP LDPC decoding
 * 
 * @return u64 
 */
unsigned ldpc_decoder::decode()
{
    u64 *vn;
    u64 *cn;

    u64 vw;
    u64 cw;

    //u64 nnz = mLdpcCode->nnz();

    //initialize
    for (u64 i = 0; i < mLdpcCode->nnz(); ++i)
    {
        mLv2c[i] = mLLRIn[mLdpcCode->c()[i]];
        mLc2v[i] = 0.0;
    }

    unsigned I = 0;
    while (I < mMaxIter)
    {
        // CN processing
        for (u64 i = 0; i < mLdpcCode->mc(); ++i)
        {
            cw = mLdpcCode->cn()[i].size();
            cn = const_cast<u64 *>(mLdpcCode->cn()[i].data());

            double tmp = 1;
            for (u64 j = 0; j < cw; ++j)
            {
                mExMsgCN[j] = 1 - 2 / (exp(mLv2c[cn[j]]) + 1); //tanh(mLv2c[cn[j]]);
                tmp *= mExMsgCN[j];
            }

            for (u64 j = 0; j < cw; ++j)
            {
                mLc2v[cn[j]] = log((mExMsgCN[j] + tmp) / (mExMsgCN[j] - tmp)); //2*atanh(tmp/mExMsgCN[j]);
            }
        }

        // VN processing and app calc
        for (u64 i = 0; i < mLdpcCode->nc(); ++i) // only transmitted bits
        {
            mLLROut[i] = mLLRIn[i];
            vw = mLdpcCode->vn()[i].size();                            // degree of VN
            vn = const_cast<u64 *>(mLdpcCode->vn()[i].data()); //neighbours of VN
            while (vw--)
            {
                mLLROut[i] += mLc2v[*vn++];
            }

            mCO[i] = (mLLROut[i] <= 0); // approx decision on ith bits

            vw = mLdpcCode->vn()[i].size();                            // degree of VN
            vn = const_cast<u64 *>(mLdpcCode->vn()[i].data()); //neighbours of VN
            while (vw--)
            {
                mLv2c[*vn] = mLLROut[i] - mLc2v[*vn];
                ++vn;
            }
        }

        if (mEarlyTerm)
        {
            if (is_codeword_legacy())
            {
                break;
            }
        }

        ++I;
    }

    return I;
}

bool ldpc_decoder::is_codeword_legacy()
{
    //calc syndrome
    bits_t s;
    for (u64 i = 0; i < mLdpcCode->mc(); i++)
    {
        s = 0;
        for (u64 j = 0; j < mLdpcCode->cn()[i].size(); j++)
        {
            s ^= mCO[mLdpcCode->c()[mLdpcCode->cn()[i][j]]];
        }

        if (s)
        {
            return false;
        }
    }

    return true;
}
} // namespace ldpc
