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
        return 0;
    }

    int ldpc_decoder_bec::decode(const vec_bits_t &channelInput)
    {
        auto &edges = mLdpcCode->H().nz_entry();

        //initialize
        for (int i = 0; i < mLdpcCode->nnz(); ++i)
        {
            mLv2c[i] = mLLRIn[edges[i].colIndex];
        }

        u32 I = 0;
        while (I < mDecoderParam.iterations)
        {
            // CN update
            for (int i = 0; i < mLdpcCode->mc(); ++i)
            {
                auto cw = mLdpcCode->H().row_neighbor()[i].size();
                auto &cn = mLdpcCode->H().row_neighbor()[i];

                mExMsgF[0] = mLv2c[cn[0].edgeIndex];
                mExMsgB[cw - 1] = mLv2c[cn[cw - 1].edgeIndex];
                for (u64 j = 1; j < cw; ++j)
                {
                    mExMsgF[j] = cn_update(mExMsgF[j - 1], mLv2c[cn[j].edgeIndex]);
                    mExMsgB[cw - 1 - j] = cn_update(mExMsgB[cw - j], mLv2c[cn[cw - j - 1].edgeIndex]);
                }

                mLc2v[cn[0].edgeIndex] = mExMsgB[1];
                mLc2v[cn[cw - 1].edgeIndex] = mExMsgF[cw - 2];
                for (u64 j = 1; j < cw - 1; ++j)
                {
                    mLc2v[cn[j].edgeIndex] = cn_update(mExMsgF[j - 1], mExMsgB[j + 1]);
                }
            }

            // VN update
            for (int i = 0; i < mLdpcCode->nc(); ++i)
            {
                // id channel output is no erasure
                // propagate output
                if (mLLRIn[i] != ERASURE)
                {
                    auto &vn = mLdpcCode->H().col_neighbor()[i]; //neighbours of VN
                    for (const auto &hi : vn)
                    {
                        mLv2c[hi.edgeIndex] = channelInput[i].value;
                    }

                    mLLROut[i] = channelInput[i].value;
                    mCO[i] = channelInput[i];
                }
                else // channel output is erasure
                {
                    auto vw = mLdpcCode->H().col_neighbor()[i].size();
                    auto &vn = mLdpcCode->H().col_neighbor()[i];

                    mExMsgF[0] = mLc2v[vn[0].edgeIndex];
                    mExMsgB[vw - 1] = mLc2v[vn[vw - 1].edgeIndex];
                    for (u64 j = 1; j < vw; ++j)
                    {
                        mExMsgF[j] = vn_update(mExMsgF[j - 1], mLc2v[vn[j].edgeIndex], channelInput[i]);
                        mExMsgB[vw - 1 - j] = vn_update(mExMsgB[vw - j], mLc2v[vn[vw - j - 1].edgeIndex], channelInput[i]);
                    }

                    mLv2c[vn[0].edgeIndex] = mExMsgB[1];
                    mLv2c[vn[vw - 1].edgeIndex] = mExMsgF[vw - 2];
                    for (u64 j = 1; j < vw - 1; ++j)
                    {
                        mLv2c[vn[j].edgeIndex] = vn_update(mExMsgF[j - 1], mExMsgB[j + 1], channelInput[i]);
                    }

                    // final decision
                    mLLROut[i] = mExMsgF[vw - 1]; //mExMsgB[0]
                    // if all incoming messages are erasures set the wrong bit
                    mCO[i] = (mLLROut[i] == ERASURE) ? -channelInput[i] : channelInput[i];
                }
            }

            if (mDecoderParam.earlyTerm)
            {
                // stop decoding when no erasures are left
                bool erasure_found = false;
                for (auto llr : mLLROut)
                {
                    if (llr == ERASURE)
                    {
                        erasure_found = true;
                        break;
                    }
                }

                if (!erasure_found)
                {
                    break;
                }
            }

            ++I;
        }

        return I;
    }
} // namespace ldpc
