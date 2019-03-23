#include "ldpc/ldpc.cuh"
#include "simulation.cuh"

#include <chrono>

using namespace ldpc;
using namespace std;

__global__ void decodeKernel(Ldpc_Decoder_cl* dec, double* llr_in, double* llr_out, const size_t& MaxIter, const bool& early_termination, size_t I)
{
	size_t* a = new size_t[10];
	printf("test :: %lu\n", I);
	/*
	while (I < MaxIter)
	{
	//decode_lyr_nodeupdate_global(llr_in);
	//decode_lyr_sumllr_global();
	//decode_lyr_appcalc_global(llr_in, llr_out);


	++I;

	if (early_termination)
	{
	if (is_codeword_global(c_out))
	break;
	}

	}
	*/
}

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

	Ldpc_Decoder_cl* cuda_dec;
	cudaMalloc(&cuda_dec, sizeof(Ldpc_Decoder_cl));
	cudaMemcpy(cuda_dec, dec, sizeof(Ldpc_Decoder_cl), cudaMemcpyHostToDevice);
	size_t iters = 0;
	decodeKernel<<<1, 1>>>(cuda_dec, llrin, llrout, 50, false, iters);

	cudaFree(cuda_dec);
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
