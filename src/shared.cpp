#include "sim/ldpcsim.h"

extern "C"
{
    void simulate(const char *codeFile, char *simFile, unsigned numThreads, bool *stopFlag, ldpc::u64 seed, ldpc::sim_results_t *res)
    {
        //setup LDPC code
        std::shared_ptr<ldpc::ldpc_code> code(new ldpc::ldpc_code(codeFile));

        // TODO - decoder parameters
        ldpc::param_map decoderParams;
        //decoderParams["iterations"] = 50;
        //decoderParams["type"] = "BP";
        //decoderParams["early_termination"] = true;

        // TODO - channel parameters
        ldpc::param_map channelParams;
        //channelParams["type"] = "AWGN";
        channelParams["seed"] = seed;
        //channelParams["x_range"] = snr;

        // TODO - simulation parameters
        ldpc::param_map simulationParams;
        simulationParams["threads"] = numThreads;
        //simulationParams["fec"] = 50;
        //simulationParams["max_frames"] = 1000000000;

        //setup sim
        ldpc::ldpc_sim sim(
            code,
            decoderParams,
            channelParams,
            simulationParams
        );

        sim.start(stopFlag);
    }

    void allocate_results(ldpc::sim_results_t *res, ldpc::u64 len)
    {
        res->fer = new double[len]();
        res->ber = new double[len]();
        res->avg_iter = new double[len]();
        res->time = new double[len]();
        res->fec = new ldpc::u64[len]();
        res->frames = new ldpc::u64[len]();
    }

    void free_results(ldpc::sim_results_t *res)
    {
        delete[] res->fer;
        delete[] res->ber;
        delete[] res->avg_iter;
        delete[] res->time;
        delete[] res->fec;
        delete[] res->frames;
    }

    ldpc::u64 calculate_rank(const char *codeFile)
    {
        ldpc::ldpc_code code(codeFile);
        return code.calc_rank();
    }
}
