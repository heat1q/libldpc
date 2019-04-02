#include "ldpc/ldpc.cuh"
#include "simulation.cuh"

using namespace ldpc;
using namespace std;


//nvcc -std=c++11 sim_cuda.cu simulation.cu ldpc/ldpc.cu ldpc/decoder.cu -o sim_cuda -arch sm_35 -rdc=true -O3
int main()
{
	Ldpc_Code_cl* code = new Ldpc_Code_cl("../src/code/test_code/code_rand_proto_3x6_400_4.txt", "../src/code/test_code/layer_rand_proto_3x6_400_4.txt", true);

	code->print_ldpc_code();

	delete code;
/*
	//set up code class on unified memory
	Ldpc_Code_cl* code_managed;
	cudaMallocManaged(&code_managed, sizeof(Ldpc_Code_cl));
	*code_managed = Ldpc_Code_cl();
	code_managed->setup_code_managed("../src/code/test_code/code_rand_proto_3x6_400_4.txt", "../src/code/test_code/layer_rand_proto_3x6_400_4.txt");
	//code_managed->setup_code_managed("../src/code/test_code/code10K.txt", "../src/code/test_code/layer10K.txt");

	//set up simulation
	//Sim_AWGN_cl sim = Sim_AWGN_cl(code_managed, "../src/sim.txt", "../src/code/test_code/map10K.txt");

	//set up decoder on unified memory
	Ldpc_Decoder_cl* dec_ufd;
	cudaMallocManaged(&dec_ufd, sizeof(Ldpc_Decoder_cl));
	dec_ufd->setup_decoder_managed(code_managed, 50, false);

	for (size_t i=0; i<code_managed->nc(); ++i)
	{
		dec_ufd->llr_in[i] = Sim_AWGN_cl::randn();
		dec_ufd->llr_out[i] = 0.0;
	}

	dec_ufd->prefetch_gpu();
	cudakernel::decode_layered<<<1,1>>>(dec_ufd);

	//TIME_PROF("CPU", dec_ufd->decode_layered_legacy(llrin, llrout, 50, false), "ms");


	//destroy decoder
	dec_ufd->destroy_dec_managed();
	cudaFree(dec_ufd);

	//destroy code
	code_managed->destroy_ldpc_code_managed();
	cudaFree(code_managed);
*/
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
