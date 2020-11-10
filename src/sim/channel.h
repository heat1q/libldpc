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

        virtual void encode_and_map();
        virtual void simulate();
        virtual void calculate_llrs();

        // Current transmitted codeword
        const vec_bits_t &codeword() const { return mCodeWord; }

        // Current encoded information word
        const vec_bits_t &infoword() const { return mInfoWord; }

    protected:
        // ptr to const ldpc_code for parameters
        std::shared_ptr<ldpc_code> mLdpcCode;

        // Holds the unique decoder for this channel
        std::shared_ptr<ldpc_decoder> mLdpcDecoder;

        // RNG engine
        std::mt19937_64 mRNG;

        // RNG for encoding information word
        std::function<bool()> mRandInfoWord;

        // Information word
        vec_bits_t mInfoWord;
        vec_bits_t mCodeWord;
    };

    class channel_awgn : public channel
    {
    public:
        channel_awgn() = default;
        channel_awgn(const std::shared_ptr<ldpc_code> &code, std::shared_ptr<ldpc_decoder> decoder, const u64 seed, const double snr);
        virtual ~channel_awgn() = default;

        /**
        * @brief Update the noise variance and reinitialize the RNG.
        * 
        * @param channelParam SNR defined as E_s / sigma^2, E_s = 1.
        */
        void set_channel_param(const double channelParam) override;

        /**
         * @brief Encodes an randomly chosen information word and maps
         * it to the channel input.
         * 
         */
        void encode_and_map() override;

        /**
        * @brief Simulate the binary-input (biAWGN) channel with noise variance sigma^2.
        * 
        */
        void simulate() override;

        /**
        * @brief Calculate the LLR values for the biAWGN channel.
        * 
        */
        void calculate_llrs() override;

    private:
        //channel i/o
        vec_double_t mX;
        vec_double_t mY;

        // signal to noise ratio, defined as 10log10(1/sigma2)
        double mSNR;
        double mSigma2;

        std::function<double()> mRandNormal;
    };

    class channel_bsc : public channel
    {
    public:
        channel_bsc() = default;
        channel_bsc(const std::shared_ptr<ldpc_code> &code, std::shared_ptr<ldpc_decoder> decoder, const u64 seed, const double epsilon);
        virtual ~channel_bsc() = default;

        /**
        * @brief Update the crossover probability and reinitialize the RNG.
        * 
        * @param channelParam espilon < 0.5
        */
        void set_channel_param(const double channelParam) override;

        /**
         * @brief Encodes an randomly chosen information word and maps
         * it to the channel input.
         * 
         */
        void encode_and_map() override;

        /**
        * @brief Simulate the Binary Symmetric Channel (BSC) with crossover probability \epsilon.
        * 
        */
        void simulate() override;

        /**
        * @brief Calculate the LLR values for the BSC.
        * 
        */
        void calculate_llrs() override;

    private:
        //channel i/o
        vec_bits_t mX;
        vec_bits_t mY;

        // crossover probability
        double mEpsilon;

        std::function<bool()> mRandBernoulli;
    };
} // namespace ldpc
