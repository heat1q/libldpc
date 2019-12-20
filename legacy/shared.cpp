#include "sim/ldpcsim.h"

extern "C"
{
    void simulate(char *codeFile, char *simFile, int numThreads, std::uint8_t *stopFlag, std::size_t seed, pgd::sim_results_t *res)
    {
        //setup LDPC code
        pgd::ldpc_code code(codeFile);

        //setup sim
        pgd::ldpc_sim sim(&code, simFile, "", numThreads, seed, res);
        sim.allocate_results();
        sim.print();

        //sim.print_file_header("", codeFile, simFile, "");
        sim.start(stopFlag);
        //sim.free_results();
    }

    void free_results(pgd::sim_results_t *res)
    {
        delete[] res->fer;
        delete[] res->ber;
        delete[] res->avg_iter;
        delete[] res->time;
        delete[] res->fec;
        delete[] res->frames;
    }

    std::size_t calculate_rank(char *codeFile)
    {
        pgd::ldpc_code code(codeFile);
        return code.calc_rank();
    }
}