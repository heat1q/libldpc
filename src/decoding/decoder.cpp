#include "decoder.h"

namespace ldpc
{
    /**
 * @brief Construct a new ldpc decoder::ldpc decoder object
 * 
 * @param pCode 
 * @param pI 
 * @param pEarlyTerm 
 */
    ldpc_decoder::ldpc_decoder(const std::shared_ptr<ldpc_code> &code, const unsigned iter, const bool earlyTerm)
        : mLdpcCode(code),
          mLv2c(code->nnz()), mLc2v(code->nnz()),
          mExMsgCN(code->max_dc()),
          mLLRIn(code->nc()), mLLROut(code->nc()),
          mSynd(code->mc()), mCO(code->nc()),
          mMaxIter(iter), mEarlyTerm(earlyTerm)
    {
    }

    /**
    * @brief Standard BP LDPC decoding
    * 
    * @return u64 
    */
    unsigned ldpc_decoder::decode()
    {
        //initialize
        for (u64 i = 0; i < mLdpcCode->nnz(); ++i)
        {
            mLv2c[i] = mLLRIn[mLdpcCode->c()[i]];
        }

        unsigned I = 0;
        while (I < mMaxIter)
        {
            // CN processing
            for (u64 i = 0; i < mLdpcCode->mc(); ++i)
            {
                auto cw = mLdpcCode->cn()[i].size();
                auto cn = const_cast<u64 *>(mLdpcCode->cn()[i].data());

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
                auto vw = mLdpcCode->vn()[i].size();                    // degree of VN
                auto vn = const_cast<u64 *>(mLdpcCode->vn()[i].data()); //neighbours of VN
                while (vw--)
                {
                    mLLROut[i] += mLc2v[*vn++];
                }

                mCO[i] = (mLLROut[i] <= 0); // approx decision on ith bits

                vw = mLdpcCode->vn()[i].size();                    // degree of VN
                vn = const_cast<u64 *>(mLdpcCode->vn()[i].data()); //neighbours of VN
                while (vw--)
                {
                    mLv2c[*vn] = mLLROut[i] - mLc2v[*vn];
                    ++vn;
                }
            }

            if (mEarlyTerm)
            {
                if (is_codeword())
                {
                    break;
                }
            }

            ++I;
        }

        return I;
    }

    bool ldpc_decoder::is_codeword()
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
