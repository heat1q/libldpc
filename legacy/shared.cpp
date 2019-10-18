#define LIB_SHARED

#include "sim/ldpcsim.h"

extern "C"

void simulate(char *codeFile, char *simFile, int numThreads)
{
    //setup LDPC code
    ldpc::ldpc_code code(codeFile);

    //setup sim
    ldpc::ldpc_sim sim(&code, simFile, "", numThreads);
    sim.print();

    //sim.print_file_header("", codeFile, simFile, "");
    sim.start();
}