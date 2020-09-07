#include "sim/ldpcsim.h"

extern "C"
{
    void simulate(char *codeFile, char *simFile, unsigned numThreads, bool *stopFlag, ldpc::u64 seed, ldpc::sim_results_t *res)
    {
        //setup LDPC code
        ldpc::ldpc_code code(codeFile);

        //setup sim
        ldpc::ldpc_sim sim(&code, simFile, "", numThreads, seed, res);
        sim.print();

        //sim.print_file_header("", codeFile, simFile, "");
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

    ldpc::u64 calculate_rank(char *codeFile)
    {
        ldpc::ldpc_code code(codeFile);
        return code.calc_rank();
    }
}
