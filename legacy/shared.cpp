#include "sim/ldpcsim.h"

extern "C"
{
    void simulate(char *codeFile, char *simFile, int numThreads, std::uint8_t *stopFlag, std::size_t seed, pgd::sim_results_t *res)
    {
        //setup LDPC code
        pgd::ldpc_code code(codeFile);

        //setup sim
        pgd::ldpc_sim sim(&code, simFile, "", numThreads, seed, res);
        sim.print();

        //sim.print_file_header("", codeFile, simFile, "");
        sim.start(stopFlag);
    }

    void allocate_results(pgd::sim_results_t *res, std::size_t len)
    {
        res->fer = new double[len]();
        res->ber = new double[len]();
        res->avg_iter = new double[len]();
        res->time = new double[len]();
        res->fec = new std::size_t[len]();
        res->frames = new std::size_t[len]();
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