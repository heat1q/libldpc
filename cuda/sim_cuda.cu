#include "ldpc/ldpc.cuh"
#include "simulation.cuh"

using namespace ldpc;
using namespace std;


//nvcc -std=c++11 sim_cuda.cu simulation.cu ldpc/ldpc.cu ldpc/decoder.cu -o sim_cuda -arch sm_35 -rdc=true -O3 -w
int main()
{
	//set up code class on unified memory
	Ldpc_Code_cl* code_dev = new Ldpc_Code_cl("../src/code/test_code/code_rand_proto_3x6_400_4.txt", "../src/code/test_code/layer_rand_proto_3x6_400_4.txt", true);

	//set up simulation
	//Sim_AWGN_cl sim = Sim_AWGN_cl(code_managed, "../src/sim.txt", "../src/code/test_code/map10K.txt");

	//set up decoder on unified memory
	Ldpc_Decoder_cl* dec_dev = new Ldpc_Decoder_cl(code_dev, 50, false, true);


	cudakernel::sim::sim_test<<<1,1>>>(dec_dev);

	cudaDeviceSynchronize();


	TIME_PROF("CPU", dec_dev->decode_layered_legacy(), "ms");
	
	TIME_PROF("CPU", dec_dev->decode_legacy(), "ms");

	//destroy decoder
	delete dec_dev;

	//destroy code
	delete code_dev;

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
