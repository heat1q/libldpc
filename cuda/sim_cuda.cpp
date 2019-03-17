#include "ldpc/ldpc.h"
#include "simulation.h"

using namespace ldpc;

int main()
{
    Ldpc_Code_cl* code = new Ldpc_Code_cl("../src/code/test_code/code_rand_proto_3x6_400_4.txt");
    code->print_ldpc_code();

    Sim_AWGN_cl* sim = new Sim_AWGN_cl(code, "../src/sim.txt", "../src/code/test_code/mapping_rand_proto_3x6_400_4.txt");
    sim->print_sim();


    delete sim;
    delete code;

    return 0;
}
