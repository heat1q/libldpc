#include "ldpc/ldpc.h"
#include "simulation.h"

#include <chrono>

using namespace ldpc;
using namespace std;

int main()
{
    /*
    Ldpc_Code_cl* code = new Ldpc_Code_cl("../src/code/test_code/code_rand_proto_3x6_400_4.txt", "../src/code/test_code/layer_rand_proto_3x6_400_4.txt");
    code->print_ldpc_code();

    Ldpc_Decoder_cl dec = Ldpc_Decoder_cl(code);
    //printVector<size_t>(dec.cn()[0], 1);


    Sim_AWGN_cl* sim = new Sim_AWGN_cl(code, "../src/sim.txt", "../src/code/test_code/mapping_rand_proto_3x6_400_4.txt");
    sim->print_sim();

    sim->start_sim();

    delete sim;


    delete code;
    */
    /*
    auto start = chrono::high_resolution_clock::now();
    for (size_t i=0; i<1000000;++i) {}
    auto elapsed = chrono::high_resolution_clock::now() - start;
    unsigned long nanoseconds = chrono::duration_cast<chrono::nanoseconds>(elapsed).count();
    cout << "Time: " << static_cast<float>(nanoseconds) << "ns" << endl;
    */

    return 0;
}
