#pragma once

#include "../decoding/decoder.h"

namespace ldpc
{
    class channel
    {
    public:
        channel() = default;
        channel(const std::shared_ptr<ldpc_code> &code, std::shared_ptr<ldpc_decoder> decoder, const u64 seed);
        virtual ~channel() = default;

        // set param and update rng stream
        virtual void set_channel_param(const double channelParam);

        virtual void simulate();
        virtual void calculate_llrs();
    protected:
        // ptr to const ldpc_code for parameters
        std::shared_ptr<ldpc_code> mLdpcCode;

        // hold the unique decoder for this channel
        std::shared_ptr<ldpc_decoder> mLdpcDecoder;

        // rng params
        const u64 mRNGSeed;
        std::mt19937_64 mRNGEngine;
    };

    class channel_awgn : public channel
    {
    public:
        channel_awgn() = default;
        channel_awgn(const std::shared_ptr<ldpc_code> &code, std::shared_ptr<ldpc_decoder> decoder, const u64 seed, const double snr);
        virtual ~channel_awgn() = default;

        void set_channel_param(const double channelParam) override;

        void simulate() override;
        void calculate_llrs() override;
    private:
        //channel i/o
        vec_double_t mX;
        vec_double_t mY;

        // signal to noise ratio, defined as 10log10(1/sigma2)
        double mSNR;
        double mSigma2;

        std::normal_distribution<double> mRandNormal;
    };

    class channel_bsc : public channel
    {
    public:
        channel_bsc() = default;
        channel_bsc(const std::shared_ptr<ldpc_code> &code, std::shared_ptr<ldpc_decoder> decoder, const u64 seed, const double epsilon);
        virtual ~channel_bsc() = default;

        void set_channel_param(const double channelParam) override;

        void simulate() override;
        void calculate_llrs() override;
    private:
        //channel i/o
        vec_bits_t mX;
        vec_bits_t mY;

        // crossover probability
        double mEpsilon;

        std::bernoulli_distribution mRandBernoulli;
    };
} // namespace ldpc
