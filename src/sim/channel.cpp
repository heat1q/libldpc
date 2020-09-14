#include "channel.h"

namespace ldpc
{
    channel::channel(const std::shared_ptr<ldpc_code> &code, std::shared_ptr<ldpc_decoder> decoder, const u64 seed)
        : mLdpcCode(code),
          mLdpcDecoder(decoder),
          mRNGSeed(seed),
          mRNGEngine(std::mt19937_64(seed))
    {
    }

    void channel::set_channel_param(const double channelParam) {}
    void channel::simulate() {}
    void channel::calculate_llrs() {}

    channel_awgn::channel_awgn(const std::shared_ptr<ldpc_code> &code, std::shared_ptr<ldpc_decoder> decoder, const u64 seed, const double snr)
        : channel(code, decoder, seed),
          mX(vec_double_t(code->nct(), 1.)), // initialize to all one, i.e. all zero cw
          mY(vec_double_t(code->nct())),
          mSNR(snr),
          mSigma2(pow(10, -snr / 10)),
          mRandNormal(std::normal_distribution<double>(0., sqrt(mSigma2)))
    {
    }

    /**
     * @brief Update the noise variance and reinitialize the RNG.
     * 
     * @param channelParam SNR defined as E_s / \sigma^2, E_s = 1.
     */
    void channel_awgn::set_channel_param(const double channelParam)
    {
        mSNR = channelParam;
        mSigma2 = pow(10, -mSNR / 10);
        mRandNormal = std::normal_distribution<double>(0., sqrt(mSigma2));
    }

    /**
     * @brief Simulate the binary-input (biAWGN) channel with noise variance \sigma^2.
     * 
     */
    void channel_awgn::simulate()
    {
        for (u64 i = 0; i < mLdpcCode->nct(); ++i)
        {
            mY[i] = mRandNormal(mRNGEngine) + mX[i]; // x \in {-1,+1}
        }
    }

    /**
     * @brief Calculate the LLR values for the biAWGN channel.
     * 
     */
    void channel_awgn::calculate_llrs()
    {
        //puncturing & shortening
        if (mLdpcCode->puncture().size() != 0)
        {
            for (auto p : mLdpcCode->puncture())
            {
                mLdpcDecoder->mLLRIn[p] = 1e-9 * (1 - 2 * (rand() % 2)); // set random signed low values, i.e. simulate erasure LLR=0
            }
        }
        if (mLdpcCode->shorten().size() != 0)
        {
            for (auto s : mLdpcCode->shorten())
            {
                mLdpcDecoder->mLLRIn[s] = 1e9; // simulate certain bit
            }
        }

        //bpsk
        for (u64 i = 0; i < mLdpcCode->nct(); ++i)
        {
            mLdpcDecoder->mLLRIn[mLdpcCode->bit_pos()[i]] = 2 * mY[i] / mSigma2;
        }
    }

    channel_bsc::channel_bsc(const std::shared_ptr<ldpc_code> &code, std::shared_ptr<ldpc_decoder> decoder, const u64 seed, const double epsilon)
        : channel(code, decoder, seed),
          mX(vec_bits_t(code->nct(), 0)), // initialize to all zero cw
          mY(vec_bits_t(code->nct())),
          mEpsilon(epsilon),
          mRandBernoulli(std::bernoulli_distribution(epsilon))
    {
    }

    /**
     * @brief Update the crossover probability and reinitialize the RNG.
     * 
     * @param channelParam \espilon < 0.5
     */
    void channel_bsc::set_channel_param(const double channelParam)
    {
        mEpsilon = channelParam;
        mRandBernoulli = std::bernoulli_distribution(mEpsilon);
    }

    /**
     * @brief Simulate the Binary Symmetric Channel (BSC) with crossover probability \epsilon.
     * 
     */
    void channel_bsc::simulate()
    {
        for (u64 i = 0; i < mLdpcCode->nct(); ++i)
        {
            mY[i] = static_cast<bits_t>(mRandBernoulli(mRNGEngine)); // generate a 1 with probability epsilon
        }
    }

    /**
     * @brief Calculate the LLR values for the BSC.
     * 
     */
    void channel_bsc::calculate_llrs()
    {
        const double delta = log((1 - mEpsilon) / mEpsilon);

        //puncturing & shortening
        if (mLdpcCode->puncture().size() != 0)
        {
            for (auto p : mLdpcCode->puncture())
            {
                mLdpcDecoder->mLLRIn[p] = delta;// * (1 - 2 * (rand() % 2)); // set random signed low values, i.e. simulate erasure LLR=0
            }
        }
        if (mLdpcCode->shorten().size() != 0)
        {
            for (auto s : mLdpcCode->shorten())
            {
                mLdpcDecoder->mLLRIn[s] = delta; // simulate certain bit
            }
        }

        for (u64 i = 0; i < mLdpcCode->nct(); ++i)
        {
            // +delta for 0; -delta for 1
            mLdpcDecoder->mLLRIn[mLdpcCode->bit_pos()[i]] = delta * (1 - 2 * mY[i]);
        }
    }

} // namespace ldpc
