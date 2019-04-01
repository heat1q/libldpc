#include "ldpc/ldpc.cuh"
#include "simulation.cuh"

using namespace ldpc;
using namespace std;


//nvcc -std=c++11 sim_cuda.cu simulation.cu ldpc/ldpc.cu ldpc/decoder.cu -o sim_cuda -arch sm_35 -rdc=true -O3
int main()
{
	//set up code class on unified memory 
	Ldpc_Code_cl* code_managed;
	cudaMallocManaged(&code_managed, sizeof(Ldpc_Code_cl));
	*code_managed = Ldpc_Code_cl();
	code_managed->setup_code_managed("../src/code/test_code/code10K.txt", "../src/code/test_code/layer10K.txt");

	//set up simulation
	//Sim_AWGN_cl sim = Sim_AWGN_cl(code_managed, "../src/sim.txt", "../src/code/test_code/map10K.txt");

	double *llrin, *llrout;
	cudaMallocManaged(&llrin, code_managed->nc()*sizeof(double));
	cudaMallocManaged(&llrout, code_managed->nc()*sizeof(double));
	for (size_t i=0; i<code_managed->nc(); ++i)
	{
		llrin[i] = Sim_AWGN_cl::randn();
		llrout[i] = 0.0;
	}

	//set up decoder on unified memory
	Ldpc_Decoder_cl* dec_ufd;
	cudaMallocManaged(&dec_ufd, sizeof(Ldpc_Decoder_cl));
	dec_ufd->setup_decoder_managed();


	TIME_PROF("GPU", dec.decode_layered(llrin, llrout, 50, false), "ms");

	TIME_PROF("CPU", dec.decode_layered_legacy(llrin, llrout, 50, false), "ms");


	cudaFree(llrin);
	cudaFree(llrout);
	
	//destroy decoder
	dec_ufd->destroy_dec_managed();
	cudaFree(dec_ufd);
	
	//destroy code
	code_managed->destroy_ldpc_code_managed();
	cudaFree(code_managed);

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