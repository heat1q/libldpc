#include "decoder.h"

namespace ldpc
{
    ldpc_decoder::ldpc_decoder(const std::shared_ptr<ldpc_code> &code,
                               const decoder_param &decoderParam)
        : ldpc_decoder_base<double>(code, decoderParam)
    {
    }

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
                auto &vn = mLdpcCode->H().col_neighbor()[i]; //neighbours of VN

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

    ldpc_decoder_bec::ldpc_decoder_bec(const std::shared_ptr<ldpc_code> &code,
                                       const decoder_param &decoderParam)
        : ldpc_decoder_base<u8>(code, decoderParam)
    {
    }

    int ldpc_decoder_bec::decode()
    {

    }
} // namespace ldpc
