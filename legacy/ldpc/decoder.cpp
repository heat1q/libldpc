#include "../sim/ldpcsim.h"

namespace ldpc
{
/*
*	Decoder
*/
//Init constructor
ldpc_decoder::ldpc_decoder(ldpc_code *pCode, std::size_t pI, bool pEarlyTerm)
    : mLdpcCode(pCode), mMaxIter(pI), mEarlyTerm(pEarlyTerm),
      mLv2c(pCode->nnz()), mLc2v(pCode->nnz()), 
      mExMsgCN(pCode->max_dc()),
      mLLRIn(pCode->nc()), mLLROut(pCode->nc()),
      mSynd(pCode->mc()), mCO(pCode->nc()),
      mIter(0), mIsCW(false)
{
}


//layered cpu decoder
std::size_t ldpc_decoder::decode_layered()
{
    std::size_t* vn;
    std::size_t* cn;

    std::size_t vw;
    std::size_t cw;

    std::size_t nnz = mLdpcCode->nnz();

    //initialize
    for (std::size_t i = 0; i < nnz; ++i)
    {
        mLv2c[i] = mLLRIn[mLdpcCode->c()[i]];
        mLc2v[i] = 0.0;
    }

    std::size_t I = 0;
    while (I < mMaxIter)
    {
        for (std::size_t l = 0; l < mLdpcCode->layers().size(); ++l)
        {
            // CN processing
            for (std::size_t i = 0; i < mLdpcCode->layers()[l].size(); ++i)
            {
                cw = mLdpcCode->cn()[mLdpcCode->layers()[l][i]].size();
                cn = const_cast<std::size_t *>(mLdpcCode->cn()[mLdpcCode->layers()[l][i]].data());
                
                double tmp = 1;
                for (std::size_t j = 0; j < cw; ++j)
                {
                    mExMsgCN[j] = 1 - 2/(exp(mLv2c[cn[j]])+1); //tanh(mLv2c[cn[j]]);
                    tmp *= mExMsgCN[j];
                }
                
                for (std::size_t j = 0; j < cw; ++j)
                {
                    mLc2v[cn[j]] = log((mExMsgCN[j]+tmp)/(mExMsgCN[j]-tmp));//2*atanh(tmp/mExMsgCN[j]);
                }
            }

            // VN processing and app calc
            for(std::size_t i = 0; i < mLdpcCode->nc(); ++i)
            {
                double tmp = mLLRIn[i];
                vw = mLdpcCode->vn()[i].size(); // degree of VN
                vn = const_cast<std::size_t *>(mLdpcCode->vn()[i].data()); //neighbours of VN
                while(vw--)
                    tmp += mLc2v[*vn++];

                mCO[i] = (mLLROut[i] <= 0); // approx decision on ith bits
                mLLROut[i] = tmp;

                vw = mLdpcCode->vn()[i].size(); // degree of VN
                vn = const_cast<std::size_t *>(mLdpcCode->vn()[i].data()); //neighbours of VN
                while(vw--)
                {
                    mLv2c[*vn] = tmp - mLc2v[*vn];
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