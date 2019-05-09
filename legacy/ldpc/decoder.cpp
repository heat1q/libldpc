#include "../sim/ldpcsim.h"

namespace ldpc
{
/*
*	Decoder device
*/
//Init constructor
ldpc_decoder::ldpc_decoder(ldpc_code *pCode, std::size_t pI, bool pEarlyTerm)
    : mLdpcCode(pCode), mMaxIter(pI), mEarlyTerm(pEarlyTerm),
      mLv2c(pCode->nnz()), mLc2v(pCode->nnz()),
      mF(pCode->max_dc()), mB(pCode->max_dc()),
      mLLRIn(pCode->nc()), mLLROut(pCode->nc()),
      mSynd(pCode->mc()), mCO(pCode->nc()),
      mIter(0), mIsCW(false)
{
}

//legacy cpu decoder
std::size_t ldpc_decoder::decode_legacy()
{
    std::size_t it;

    std::size_t *vn;
    std::size_t *cn;

    std::size_t vw;
    std::size_t cw;

       
    //initialize with llrs
    for (std::size_t i = 0; i < mLdpcCode->nnz(); i++)
    {
        mLv2c[i] = mLLRIn[mLdpcCode->c()[i]];
    }

    it = 0;
    while (it < mMaxIter)
    {
        for (std::size_t i = 0; i < mLdpcCode->mc(); i++)
        {
            cw = mLdpcCode->cn()[i].size();
            cn = const_cast<std::size_t *>(mLdpcCode->cn()[i].data());
            mF[0] = mLv2c[*cn];
            mB[cw - 1] = mLv2c[*(cn + cw - 1)];
            for (std::size_t j = 1; j < cw; j++)
            {
                mF[j] = jacobian(mF[j - 1], mLv2c[*(cn + j)]);
                mB[cw - 1 - j] = jacobian(mB[cw - j], mLv2c[*(cn + cw - j - 1)]);
            }

            mLc2v[*cn] = mB[1];
            mLc2v[*(cn + cw - 1)] = mF[cw - 2];
            for (std::size_t j = 1; j < cw - 1; j++)
            {
                mLc2v[*(cn + j)] = jacobian(mF[j - 1], mB[j + 1]);
            }
        }

        // VN node processing
        for (std::size_t i = 0; i < mLdpcCode->nc(); i++)
        {
            double tmp = mLLRIn[i];
            vw = mLdpcCode->vn()[i].size();
            vn = const_cast<std::size_t *>(mLdpcCode->vn()[i].data());
            while (vw--)
            {
                tmp += mLc2v[*vn++];
            }
            vn = const_cast<std::size_t *>(mLdpcCode->vn()[i].data());
            vw = mLdpcCode->vn()[i].size();
            while (vw--)
            {
                mLv2c[*vn] = tmp - mLc2v[*vn];
                vn++;
            }
        }

        // app calculation
        for (std::size_t i = 0; i < mLdpcCode->nc(); i++)
        {
            mLLROut[i] = mLLRIn[i];
            vn = const_cast<std::size_t *>(mLdpcCode->vn()[i].data());
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

bool ldpc_decoder::is_codeword_legacy()
{
    //calc syndrome
    bits_t s;
    for (std::size_t i = 0; i < mLdpcCode->mc(); i++)
    {
        s = 0;
        for (std::size_t j = 0; j < mLdpcCode->cw()[i]; j++)
            s ^= mCO[mLdpcCode->c()[mLdpcCode->cn()[i][j]]];

        if (s)
        {
            return false;
        }
    }

    return true;
}
} // namespace ldpc