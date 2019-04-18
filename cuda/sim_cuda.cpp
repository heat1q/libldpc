#include "ldpcsim.h"
#include "device/vectormgd.h"

using namespace ldpc;

// /usr/local/cuda-9.2/bin/nvcc -x cu -std=c++11 sim_cuda.cpp ldpcsim.cpp ldpc/ldpc.cpp ldpc/decoder.cpp device/cudamgd.cpp device/kernel.cpp device/ldpcsimdevice.cpp  -o sim_cuda -arch sm_35 -rdc=true -O3 -w
int main()
{
/*
    //set up code class on unified memory
    ldpc_code* code = new ldpc_code("../src/code/test_code/code_rand_proto_3x6_400_4.txt", "../src/code/test_code/layer_rand_proto_3x6_400_4.txt", true);

	cudamgd_ptr<ldpc_code> code(code_dev);

    //set up simulation
    ldpc_sim sim(code_dev, "../src/sim.txt", "../src/code/test_code/mapping_rand_proto_3x6_400_4.txt");
    sim.print();
    sim.start();
*/


	ldpc_code_device code_dev(
		"../src/code/test_code/code_rand_proto_3x6_400_4.txt"
		, "../src/code/test_code/layer_rand_proto_3x6_400_4.txt"
	);

	cudamgd_ptr<ldpc_sim_device> sim(
		ldpc_sim_device(
			code_dev
			, "../src/sim.txt"
			, "../src/code/test_code/mapping_rand_proto_3x6_400_4.txt"
		)
	);
	//sim->print();
	sim->start();

/*
    //set up decoder on unified memory
	try
	{
		cudamgd_ptr<ldpc_decoder_device> dec_dev(
			ldpc_decoder_device(code_dev, 50, false)
		);

		for (size_t i = 0; i < code_dev.nc(); ++i) {
			dec_dev->mLLRIn[i] = ldpc_sim::randn();
		}

		dec_dev->decode_layered();
		for (size_t i = 0; i < code_dev.nc(); ++i) {
			double a = dec_dev->mLLRIn[i];
		}
		dec_dev->decode_layered();

		//TIME_PROF("GPU Layered", dec_dev->decode_layered(), "ms");
		//TIME_PROF("CPU Layered", dec_dev->decode_layered_legacy(), "ms");
		//TIME_PROF("CPU Legacy", dec_dev->decode_legacy(), "ms");
	}
	catch(...)
	{
		std::cout << "Error ldpc_decoder_device" << '\n';
	}
    return 0;
*/
}
