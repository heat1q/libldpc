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
          mCO(code->nc()),
          mCNApprox(jacobian),
          mDecoderParam(decoderParam)
    {
        set_approximation();
    }

    void ldpc_decoder::set_approximation()
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
    int ldpc_decoder::decode()
    {
        auto &edges = mLdpcCode->H().nz_entry();

        //initialize
        for (int i = 0; i < mLdpcCode->nnz(); ++i)
        {
            mLv2c[i] = mLLRIn[edges[i].colIndex];
        }

        unsigned I = 0;
        while (I < mDecoderParam.iterations)
        {
            // CN processing
            for (int i = 0; i < mLdpcCode->mc(); ++i)
            {
                auto cw = mLdpcCode->H().row_neighbor()[i].size();
                auto &cn = mLdpcCode->H().row_neighbor()[i];

                // J. Chen et al. “Reduced-Complexity Decoding of LDPC Codes”
                mExMsgF[0] = mLv2c[cn[0].edgeIndex];
                mExMsgB[cw - 1] = mLv2c[cn[cw - 1].edgeIndex];
                for (u64 j = 1; j < cw; ++j)
                {
                    mExMsgF[j] = mCNApprox(mExMsgF[j - 1], mLv2c[cn[j].edgeIndex]);
                    mExMsgB[cw - 1 - j] = mCNApprox(mExMsgB[cw - j], mLv2c[cn[cw - j - 1].edgeIndex]);
                }

                mLc2v[cn[0].edgeIndex] = mExMsgB[1];
                mLc2v[cn[cw - 1].edgeIndex] = mExMsgF[cw - 2];
                for (u64 j = 1; j < cw - 1; ++j)
                {
                    mLc2v[cn[j].edgeIndex] = mCNApprox(mExMsgF[j - 1], mExMsgB[j + 1]);
                }
            }

            // VN processing and app calc
            for (int i = 0; i < mLdpcCode->nc(); ++i) // only transmitted bits
            {
                mLLROut[i] = mLLRIn[i];
                auto &vn = mLdpcCode->H().col_neighbor()[i];       //neighbours of VN

                for (const auto &hi : vn)
                {
                    mLLROut[i] += mLc2v[hi.edgeIndex];
                }

                mCO[i] = (mLLROut[i] <= 0); // approx decision on ith bits

                for (const auto &hi : vn)
                {
                    mLv2c[hi.edgeIndex] = mLLROut[i] - mLc2v[hi.edgeIndex];
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
        for (int i = 0; i < mLdpcCode->mc(); i++)
        {
            s = 0;
            for (const auto &hj : mLdpcCode->H().row_neighbor()[i])
            {
                s += mCO[hj.nodeIndex];
            }

            if (s != 0)
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
