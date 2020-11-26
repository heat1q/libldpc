// Shared library wrapper
#include "sim/ldpcsim.h"

static std::shared_ptr<ldpc::ldpc_code> ldpcCode;
static std::shared_ptr<ldpc::ldpc_decoder> ldpcDecoder;

using namespace ldpc;

extern "C"
{
    void ldpc_setup(const char *pcFile, 
                    const char *genFile, 
                    int *n, 
                    int *m,
                    int *nct,
                    int *mct)
    {
        ldpcCode = std::make_shared<ldpc::ldpc_code>(pcFile, genFile);
        decoder_param decoderParams;
        decoderParams.type = "";
        ldpcDecoder = std::make_shared<ldpc::ldpc_decoder>(ldpcCode, decoderParams);
        *n = ldpcCode->nc(); *m = ldpcCode->mc();
        *nct = ldpcCode->nct(); *mct = ldpcCode->mct();
    }

    void simulate(decoder_param decoderParams, channel_param channelParam, simulation_param simParam, sim_results_t *results, bool *stopFlag)
    {
        ldpc_sim sim(ldpcCode, decoderParams, channelParam, simParam, results);
        sim.start(stopFlag);
    }

    int calculate_rank()
    {
        return ldpcCode->H().rank();
    }

    void encode(uint8_t *infoWord, uint8_t *codeWord)
    {
        vec_bits_t u(infoWord, infoWord + ldpcCode->kct());
        auto cw = ldpcCode->G().multiply_left(u);
        for (int i = 0; i < ldpcCode->nct(); ++i)
        {
            codeWord[i] = cw[ldpcCode->bit_pos()[i]].value;
        }
    }

    int decode(ldpc::decoder_param decoderParams, double *llr, double *llrOut)
    {
        ldpcDecoder->set_param(decoderParams);
        ldpc::vec_double_t llrIn(ldpcCode->nc(), 0.0);
        for (int i = 0; i < ldpcCode->nct(); ++i)
        {
            llrIn[ldpcCode->bit_pos()[i]] = llr[i];
        }
        ldpcDecoder->llr_in(llrIn);


        int iter = ldpcDecoder->decode();
        for (int i = 0; i < ldpcCode->nct(); ++i)
        {
            llrOut[i] = ldpcDecoder->llr_out()[ldpcCode->bit_pos()[i]];
        }
        
        return iter;
    }

    void syndrome(uint8_t* word, uint8_t* syndrome)
    {
        vec_bits_t v(word, word + ldpcCode->nc());
        
        auto s = ldpcCode->H().multiply_right(v);

        for (u64 i = 0; i < s.size(); ++i)
        {
            syndrome[i] = s[i].value;
        }
    }
}
