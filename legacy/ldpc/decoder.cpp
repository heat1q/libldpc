#include "../sim/ldpcsim.h"

namespace ldpc
{
/*
*	Decoder device
*/
//Init constructor
ldpc_decoder::ldpc_decoder(ldpc_code *pCode, std::size_t pI, bool pEarlyTerm)
    : mLdpcCode(pCode), mMaxIter(pI), mEarlyTerm(pEarlyTerm),
      mLv2c(pCode->layers().size() * pCode->nnz()), mLc2v(pCode->layers().size() * pCode->nnz()), 
      mLc2vPre(pCode->layers().size() * pCode->nnz()), mLSum(pCode->nnz()),
      mF(pCode->layers().size() * pCode->max_dc()), mB(pCode->layers().size() * pCode->max_dc()),
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
        mLSum[i] = 0.0;
        for (std::size_t l = 0; l < mLdpcCode->nl(); ++l)
        {
            mLc2v[l*nnz+i] = 0.0;
            mLv2c[l*nnz+i] = 0.0;
            mLc2vPre[l*nnz+i] = 0.0;
        }
    }

    std::size_t I = 0;
    while (I < mMaxIter)
    {
        for (std::size_t l = 0; l < mLdpcCode->nl(); ++l)
        {
            /* VN node intialization */
            for(std::size_t i = 0; i < mLdpcCode->nc(); i++)
            {
                double tmp = mLLRIn[i];
                vw = mLdpcCode->vn()[i].size();
                vn = const_cast<std::size_t *>(mLdpcCode->vn()[i].data());
                while(vw--)
                    tmp += mLSum[*vn++];

                vn = const_cast<std::size_t *>(mLdpcCode->vn()[i].data());
                vw = mLdpcCode->vn()[i].size();
                while(vw--)
                {
                    mLv2c[l*nnz + *vn] = tmp - mLc2v[l*nnz+*vn];
                    ++vn;
                }
            }

            /* CN node processing */
            for(std::size_t i = 0; i < mLdpcCode->layers()[l].size(); i++)
            {
                cw = mLdpcCode->cn()[mLdpcCode->layers()[l][i]].size();
                cn = const_cast<std::size_t *>(mLdpcCode->cn()[mLdpcCode->layers()[l][i]].data());
                mF[l*mLdpcCode->max_dc() + 0] = mLv2c[l*nnz + *cn];
                mB[l*mLdpcCode->max_dc() + cw-1] = mLv2c[l*nnz + *(cn+cw-1)];
                for(std::size_t j = 1; j < cw; j++)
                {
                    mF[l*mLdpcCode->max_dc() + j] = jacobian(mF[l*mLdpcCode->max_dc() + j-1], mLv2c[l*nnz + *(cn+j)]);
                    mB[l*mLdpcCode->max_dc() + cw-1-j] = jacobian(mB[l*mLdpcCode->max_dc() + cw-j], mLv2c[l*nnz + *(cn + cw-j-1)]);
                }

                mLc2v[l*nnz + *cn] = mB[l*mLdpcCode->max_dc() + 1];
                mLc2v[l*nnz + *(cn+cw-1)] = mF[l*mLdpcCode->max_dc() + cw-2];

                for(std::size_t j = 1; j < cw-1; j++)
                    mLc2v[l*nnz + *(cn+j)] = jacobian(mF[l*mLdpcCode->max_dc() + j-1], mB[l*mLdpcCode->max_dc() + j+1]);
            }

            //update the llr sum of layers, by replacing old llr of lyr l with new value
            for (std::size_t i = 0; i < nnz; ++i)
            {
                mLSum[i] += mLc2v[l*nnz + i] - mLc2vPre[l*nnz + i];
                mLc2vPre[l*nnz + i] = mLc2v[l*nnz + i];
            }

            // app calculation
            for(std::size_t i = 0; i < mLdpcCode->nc(); ++i)
            {
                mLLROut[i] = mLLRIn[i];
                vn = const_cast<std::size_t *>(mLdpcCode->vn()[i].data());
                vw = mLdpcCode->vn()[i].size();
                while(vw--)
                    mLLROut[i] += mLSum[*vn++];
                mCO[i] = (mLLROut[i] <= 0);
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