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
ldpc_decoder::ldpc_decoder(ldpc_code *pCode, ldpc_sim *pSim, std::int16_t pI, bool pEarlyTerm)
    : mLdpcCode(pCode), mSim(pSim), mMaxIter(pI), mEarlyTerm(pEarlyTerm),
      mLv2c(pCode->nnz()), mLc2v(pCode->nnz()),
      mExMsgCN(pCode->max_dc()),
      mLLRIn(pCode->nc()), mLLROut(pCode->nc()),
      mSynd(pCode->mc()), mCO(pCode->nc())
{
}

void ldpc_decoder::calc_llrs(const vec_double_t &y, double sigma2)
{
    //puncturing & shortening
    if (mLdpcCode->puncture().size() != 0)
    {
        for (auto p: mLdpcCode->puncture())
        {
            mLLRIn[p] = +99999.9;
        }
    }
    if (mLdpcCode->shorten().size() != 0)
    {
        for (auto s: mLdpcCode->shorten())
        {
            mLLRIn[s] = -99999.9;
        }
    }

    //bpsk
    for (std::size_t i = 0; i < mLdpcCode->nct(); ++i)
    {
        mLLRIn[mSim->bits_pos()[i]] = 2 * y[i] / sigma2;
    }
}

/**
 * @brief Standard BP LDPC decoding
 * 
 * @return std::size_t 
 */
std::int16_t ldpc_decoder::decode()
{
    std::size_t *vn;
    std::size_t *cn;

    std::size_t vw;
    std::size_t cw;

    std::size_t nnz = mLdpcCode->nnz();

    //initialize
    for (std::size_t i = 0; i < nnz; ++i)
    {
        mLv2c[i] = mLLRIn[mLdpcCode->c()[i]];
        mLc2v[i] = 0.0;
    }

    std::int16_t I = 0;
    while (I < mMaxIter)
    {
        // CN processing
        for (std::size_t i = 0; i < mLdpcCode->mc(); ++i)
        {
            cw = mLdpcCode->cn()[i].size();
            cn = const_cast<std::size_t *>(mLdpcCode->cn()[i].data());

            double tmp = 1;
            for (std::size_t j = 0; j < cw; ++j)
            {
                mExMsgCN[j] = 1 - 2 / (exp(mLv2c[cn[j]]) + 1); //tanh(mLv2c[cn[j]]);
                tmp *= mExMsgCN[j];
            }

            for (std::size_t j = 0; j < cw; ++j)
            {
                mLc2v[cn[j]] = log((mExMsgCN[j] + tmp) / (mExMsgCN[j] - tmp)); //2*atanh(tmp/mExMsgCN[j]);
            }
        }

        // VN processing and app calc
        for (std::size_t i = 0; i < mLdpcCode->nc(); ++i)
        {
            mLLROut[i] = mLLRIn[i];
            vw = mLdpcCode->vn()[i].size();                            // degree of VN
            vn = const_cast<std::size_t *>(mLdpcCode->vn()[i].data()); //neighbours of VN
            while (vw--)
            {
                mLLROut[i] += mLc2v[*vn++];
            }

            mCO[i] = (mLLROut[i] <= 0); // approx decision on ith bits

            vw = mLdpcCode->vn()[i].size();                            // degree of VN
            vn = const_cast<std::size_t *>(mLdpcCode->vn()[i].data()); //neighbours of VN
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
    for (std::size_t i = 0; i < mLdpcCode->mc(); i++)
    {
        s = 0;
        for (std::size_t j = 0; j < mLdpcCode->cn()[i].size(); j++)
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