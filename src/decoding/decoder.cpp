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
    ldpc_decoder::ldpc_decoder(const std::shared_ptr<ldpc_code> &code, const decoder_param &decoderParam)
        : mLdpcCode(code),
          mLv2c(code->nnz()), mLc2v(code->nnz()),
          mExMsgF(code->max_dc()), mExMsgB(code->max_dc()),
          mLLRIn(code->nc()), mLLROut(code->nc()),
          mSynd(code->mc()), mCO(code->nc()),
          mCNApprox(jacobian),
          mDecoderParam(decoderParam)
    {
        if (mDecoderParam.type == std::string("BP_MS"))
        {
            mCNApprox = minsum;
        }
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
        while (I < mDecoderParam.iterations)
        {
            // CN processing
            for (u64 i = 0; i < mLdpcCode->mc(); ++i)
            {
                auto cw = mLdpcCode->cn()[i].size();
                auto &cn = mLdpcCode->cn()[i];

                // J. Chen et al. “Reduced-Complexity Decoding of LDPC Codes”
                mExMsgF[0] = mLv2c[cn[0]];
                mExMsgB[cw - 1] = mLv2c[cn[cw - 1]];
                for (u64 j = 1; j < cw; ++j)
                {
                    mExMsgF[j] = mCNApprox(mExMsgF[j - 1], mLv2c[cn[j]]);
                    mExMsgB[cw - 1 - j] = mCNApprox(mExMsgB[cw - j], mLv2c[cn[cw - j - 1]]);
                }

                mLc2v[cn[0]] = mExMsgB[1];
                mLc2v[cn[cw - 1]] = mExMsgF[cw - 2];
                for (u64 j = 1; j < cw - 1; ++j)
                {
                    mLc2v[cn[j]] = mCNApprox(mExMsgF[j - 1], mExMsgB[j + 1]);
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

            if (mDecoderParam.earlyTerm)
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

    constexpr int ldpc_decoder::sign(const double x)
    {
        return (1 - 2 * static_cast<int>(std::signbit(x)));
    }

    constexpr double ldpc_decoder::jacobian(const double x, const double y)
    {
        return sign(x) * sign(y) * std::min(std::abs(x), std::abs(y)) + std::log((1 + std::exp(-std::abs(x + y))) / (1 + std::exp(-std::abs(x - y))));
    }

    constexpr double ldpc_decoder::minsum(const double x, const double y)
    {
        return sign(x) * sign(y) * std::min(std::abs(x), std::abs(y));
    }
} // namespace ldpc
