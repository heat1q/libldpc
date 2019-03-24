#include "ldpc/ldpc.cuh"
#include "simulation.cuh"

using namespace ldpc;
using namespace std;


//nvcc -std=c++11 sim_cuda.cu simulation.cu ldpc/ldpc.cu ldpc/decoder.cu -o sim_cuda
// ??? -arch sm_35 -rdc=true
int main()
{
	Ldpc_Code_cl* code_managed;// = new Ldpc_Code_cl("../src/code/test_code/code_rand_proto_3x6_400_4.txt", "../src/code/test_code/layer_rand_proto_3x6_400_4.txt");
	cudaMallocManaged(&code_managed, sizeof(Ldpc_Code_cl));
	code_managed->setup_code_managed("../src/code/test_code/code_rand_proto_3x6_400_4.txt", "../src/code/test_code/layer_rand_proto_3x6_400_4.txt");
	printf("Host :: %lu\n", code_managed->vw()[1087]);

	Ldpc_Decoder_cl** dec_ptr;
	cudaMallocManaged(&dec_ptr, sizeof(Ldpc_Decoder_cl*));

	cudakernel::setup_decoder<<<1, 1>>>(code_managed, dec_ptr);
	cudakernel::destroy_decoder<<<1, 1>>>(dec_ptr);

	code_managed->destroy_ldpc_code_managed();
	cudaFree(code_managed);
	cudaFree(dec_ptr);
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
