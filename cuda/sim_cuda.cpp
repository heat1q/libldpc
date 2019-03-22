#include "ldpc/ldpc.h"
#include "simulation.h"

#include <chrono>

using namespace ldpc;
using namespace std;

int main()
{
    Ldpc_Code_cl* code = new Ldpc_Code_cl("../src/code/test_code/code_rand_proto_3x6_400_4.txt", "../src/code/test_code/layer_rand_proto_3x6_400_4.txt");
    //code->print_ldpc_code();

    //Sim_AWGN_cl* sim = new Sim_AWGN_cl(code, "../src/sim.txt", "../src/code/test_code/mapping_rand_proto_3x6_400_4.txt");
    //sim->print_sim();
    //sim->start_sim();

    Ldpc_Decoder_cl* dec = new Ldpc_Decoder_cl(code);
    double* llrin = new double[code->nc()]();
    double* llrout = new double[code->nc()]();
    for (size_t i=0; i<code->nc(); ++i)
        llrin[i] = Sim_AWGN_cl::randn();


    auto start = chrono::high_resolution_clock::now();
    dec->decode_layered(llrin, llrout, 50, false);
    auto elapsed = chrono::high_resolution_clock::now() - start;
    cout << "Time: " << static_cast<float>(chrono::duration_cast<chrono::microseconds>(elapsed).count()) << "us" << endl;

    start = chrono::high_resolution_clock::now();
    dec->decode(llrin, llrout, 50, false);
    elapsed = chrono::high_resolution_clock::now() - start;
    cout << "Time: " << static_cast<float>(chrono::duration_cast<chrono::microseconds>(elapsed).count()) << "us" << endl;





    delete dec;
    //delete sim;
    delete code;
    return 0;
}

//tmpl fcts need definition in each file?
template<typename T> void ldpc::printVector(T *x, const size_t &l)
{
    cout << "[";
    for (size_t i = 0; i < l-1; ++i)
        cout << x[i] << " ";
    cout << x[l-1] << "]";
}

