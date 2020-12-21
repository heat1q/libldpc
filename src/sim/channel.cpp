#include "channel.h"

namespace ldpc
{
    channel::channel(const std::shared_ptr<ldpc_code> &code,
                     const decoder_param &decoderParams,
                     const u64 seed)
        : mLdpcCode(code),
          mLdpcDecoder(std::make_shared<ldpc_decoder>(code, decoderParams)),
          mRNG(seed),
          mRandInfoWord(std::bind(std::bernoulli_distribution(0.5), std::mt19937_64(seed << 1))),
          mInfoWord(vec_bits_t(code->kc(), 0)),
          mCodeWord(vec_bits_t(code->nc(), 0))
    {
    }

    void channel::set_channel_param(const double channelParam) {}
    void channel::encode_and_map() {}
    void channel::simulate() {}
    void channel::calculate_llrs() {}

    channel_awgn::channel_awgn(const std::shared_ptr<ldpc_code> &code,
                               const decoder_param &decoderParams,
                               const u64 seed,
                               const double snr)
        : channel(code, decoderParams, seed),
          mX(vec_double_t(code->nct(), 1.)), // initialize to all one, i.e. all zero cw
          mY(vec_double_t(code->nct())),
          mSNR(snr),
          mSigma2(pow(10, -snr / 10)),
          mRandNormal(std::bind(std::normal_distribution<double>(0., sqrt(mSigma2)), mRNG))
    {
    }

    void channel_awgn::set_channel_param(const double channelParam)
    {
        mSNR = channelParam;
        mSigma2 = pow(10, -mSNR / 10);
        mRandNormal = std::bind(std::normal_distribution<double>(0., sqrt(mSigma2)), mRNG);
    }

    void channel_awgn::encode_and_map()
    {
        for (auto &u : mInfoWord)
        {
            u = mRandInfoWord();
        }

        mLdpcCode->G().multiply_left(mInfoWord, mCodeWord);

        // only select transmitted codeword bits
        for (int i = 0; i < mLdpcCode->nct(); ++i)
        {
            // 0 --> +1
            // 1 --> -1
            mX[i] = 1 - (2 * mCodeWord[mLdpcCode->bit_pos()[i]].value);
        }
    }

    void channel_awgn::simulate()
    {
        for (int i = 0; i < mLdpcCode->nct(); ++i)
        {
            mY[i] = mRandNormal() + mX[i]; // x \in {-1,+1}
        }
    }

    void channel_awgn::calculate_llrs()
    {
        //puncturing & shortening
        if (mLdpcCode->puncture().size() != 0)
        {
            for (auto p : mLdpcCode->puncture())
            {
                mLdpcDecoder->mLLRIn[p] = 0.0; // punctured bits = erasure LLR=0
            }
        }
        if (mLdpcCode->shorten().size() != 0)
        {
            for (auto s : mLdpcCode->shorten())
            {
                mLdpcDecoder->mLLRIn[s] = 99999.9; // simulate certain bit
            }
        }

        //bpsk
        for (int i = 0; i < mLdpcCode->nct(); ++i)
        {
            mLdpcDecoder->mLLRIn[mLdpcCode->bit_pos()[i]] = 2 * mY[i] / mSigma2;
        }
    }

    channel_bsc::channel_bsc(const std::shared_ptr<ldpc_code> &code,
                             const decoder_param &decoderParams,
                             const u64 seed,
                             const double epsilon)
        : channel(code, decoderParams, seed),
          mX(vec_bits_t(code->nct(), 0)), // initialize to all zero cw
          mY(vec_bits_t(code->nct())),
          mEpsilon(epsilon),
          mRandBernoulli(std::bind(std::bernoulli_distribution(epsilon), mRNG))
    {
    }

    void channel_bsc::encode_and_map()
    {
        for (auto &u : mInfoWord)
        {
            u = mRandInfoWord();
        }

        mLdpcCode->G().multiply_left(mInfoWord, mCodeWord);

        // only select transmitted codeword bits
        for (int i = 0; i < mLdpcCode->nct(); ++i)
        {
            mX[i] = mCodeWord[mLdpcCode->bit_pos()[i]];
        }
    }

    void channel_bsc::set_channel_param(const double channelParam)
    {
        mEpsilon = channelParam;
        mRandBernoulli = std::bind(std::bernoulli_distribution(mEpsilon), mRNG);
    }

    void channel_bsc::simulate()
    {
        for (int i = 0; i < mLdpcCode->nct(); ++i)
        {
            mY[i] = mX[i] + mRandBernoulli(); // flip bit with probability epsilon
        }
    }

    void channel_bsc::calculate_llrs()
    {
        const double delta = log((1 - mEpsilon) / mEpsilon);

        //puncturing & shortening
        if (mLdpcCode->puncture().size() != 0)
        {
            for (auto p : mLdpcCode->puncture())
            {
                mLdpcDecoder->mLLRIn[p] = 0.0; // simulate erasure LLR=0
            }
        }
        if (mLdpcCode->shorten().size() != 0)
        {
            for (auto s : mLdpcCode->shorten())
            {
                mLdpcDecoder->mLLRIn[s] = delta; // simulate certain bit
            }
        }

        for (int i = 0; i < mLdpcCode->nct(); ++i)
        {
            // +delta for 0; -delta for 1
            mLdpcDecoder->mLLRIn[mLdpcCode->bit_pos()[i]] = delta * (1 - 2 * mY[i].value);
        }
    }

    channel_bec::channel_bec(const std::shared_ptr<ldpc_code> &code,
                             const decoder_param &decoderParams,
                             const u64 seed,
                             const double epsilon)
        : channel(code, decoderParams, seed),
          mLdpcDecoder(std::make_shared<ldpc_decoder_bec>(code, decoderParams)), // bec decoder
          mX(vec_bits_t(code->nct(), 0)), // initialize to all zero cw
          mY(std::vector<u8>(code->nct())),
          mEpsilon(epsilon),
          mRandBernoulli(std::bind(std::bernoulli_distribution(epsilon), mRNG))
    {
    }

    void channel_bec::encode_and_map()
    {
        for (auto &u : mInfoWord)
        {
            u = mRandInfoWord();
        }

        mLdpcCode->G().multiply_left(mInfoWord, mCodeWord);

        // only select transmitted codeword bits
        for (int i = 0; i < mLdpcCode->nct(); ++i)
        {
            mX[i] = mCodeWord[mLdpcCode->bit_pos()[i]];
        }
    }

    void channel_bec::set_channel_param(const double channelParam)
    {
        mEpsilon = channelParam;
        mRandBernoulli = std::bind(std::bernoulli_distribution(mEpsilon), mRNG);
    }

    void channel_bec::simulate()
    {
        for (int i = 0; i < mLdpcCode->nct(); ++i)
        {
            mY[i] = mRandBernoulli() ? ERASURE : mX[i].value; // erasure with probability epsilon
        }
    }

    void channel_bec::calculate_llrs()
    {
        //puncturing & shortening
        if (mLdpcCode->puncture().size() != 0)
        {
            for (auto p : mLdpcCode->puncture())
            {
                mLdpcDecoder->mLLRIn[p] = ERASURE; // simulate erasure LLR=0
            }
        }
        if (mLdpcCode->shorten().size() != 0)
        {
            for (auto s : mLdpcCode->shorten())
            {
                mLdpcDecoder->mLLRIn[s] = mX[s].value; // simulate certain bit
            }
        }

        for (int i = 0; i < mLdpcCode->nct(); ++i)
        {
            mLdpcDecoder->mLLRIn[mLdpcCode->bit_pos()[i]] = mY[i];
        }
    }
} // namespace ldpc
