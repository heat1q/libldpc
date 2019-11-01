#include "sim/ldpcsim.h"

extern "C"

void simulate(char *codeFile, char *simFile, int numThreads, std::uint8_t *stopFlag, std::size_t seed)
{
    //setup LDPC code
    ldpc::ldpc_code code(codeFile);

    //setup sim
    ldpc::ldpc_sim sim(&code, simFile, "", numThreads, seed);
    sim.print();

    //sim.print_file_header("", codeFile, simFile, "");
    sim.start(stopFlag);
}