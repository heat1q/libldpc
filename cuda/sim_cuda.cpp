#include "ldpcsim.h"

// /usr/local/cuda-9.2/bin/nvcc -x cu -std=c++11 sim_cuda.cpp ldpc/ldpc.cpp ldpc/decoder.cpp device/kernel.cpp device/ldpcsimdevice.cpp  -o sim_cuda -arch sm_35 -rdc=true -O3 -w
int main(int argc, char *argv[])
{
	bool abort = false;

	std::string codeFile;
	std::string simFile;
	std::string mapFile;
	std::string layerFile;
	ldpc::labels_t numThreads = 1;

	if (argc == 7)
	{
		for (int i = 0; i < (argc - 1) / 2; ++i)
		{
			if (strcmp(argv[2 * i + 1], "-code") == 0)
			{
				codeFile.assign(argv[2 * i + 2]);
			}
			else if (strcmp(argv[2 * i + 1], "-sim") == 0)
			{
				simFile.assign(argv[2 * i + 2]);
			}
			else if (strcmp(argv[2 * i + 1], "-map") == 0)
			{
				mapFile.assign(argv[2 * i + 2]);
			}
			else
			{
				abort = true;
			}
		}
	}
	else if (argc == 9 || argc == 11)
	{
		for (int i = 0; i < (argc - 1) / 2; ++i)
		{
			if (strcmp(argv[2 * i + 1], "-code") == 0)
			{
				codeFile.assign(argv[2 * i + 2]);
			}
			else if (strcmp(argv[2 * i + 1], "-sim") == 0)
			{
				simFile.assign(argv[2 * i + 2]);
			}
			else if (strcmp(argv[2 * i + 1], "-map") == 0)
			{
				mapFile.assign(argv[2 * i + 2]);
			}
			else if (strcmp(argv[2 * i + 1], "-layer") == 0)
			{
				layerFile.assign(argv[2 * i + 2]);
			}
			else if (strcmp(argv[2 * i + 1], "-threads") == 0)
			{
				numThreads = static_cast<unsigned short>(atoi(argv[2 * i + 2]));
			}
			else
			{
				abort = true;
			}
		}
	}
	else
	{
		abort = true;
	}

	std::cout << "==================================== LDPC Simulation ===================================\n";
	for (std::size_t i = 0; i < argc; ++i)
	{
		std::cout << argv[i] << " ";
	}
	std::cout << "\n\n DFLAGS: SIM_NUM_BITS=" << SIM_NUM_BITS << ", DEC_MAX_DC=" << DEC_MAX_DC
			  << ", NUMK_THREADS=" << NUMK_THREADS;
#ifdef LOG_FRAME_TIME
	std::cout << ", LOG_FRAME_TIME";
#endif
#ifdef SHORT_LOG
	std::cout << ", SHORT_LOG";
#endif
	std::cout << "\n";
	if (abort)
	{
		std::cout << "==================================== USAGE REMINDER ====================================\n";
		std::cout << "                    ./Main  -code CodeFile -sim SimFile -map MappingFile                 \n";
		std::cout << "                  optional: -layer LayerFile                                             \n";
		std::cout << "                            -threads Num                                                 \n";
		std::cout << "                                                                                         \n";
		std::cout << "                  CodeFile: Name of the code file                                        \n";
		std::cout << "               CodeMapFile: Name of mapping file                                         \n";
		std::cout << "                   SimFile: Name of simulation file                                      \n";
		std::cout << "                                                                                         \n";
		std::cout << "                 LayerFile: Name of code layer file for layered decoding                 \n";
		std::cout << "                       Num: Number of frames processed in parallel (check GPU Load)      \n";
		std::cout << "                            (default: 1)                                                 \n";
		std::cout << "========================================================================================" << std::endl;
		exit(EXIT_FAILURE);
	}

	ldpc::cuda_ptr<ldpc::ldpc_code_device> code_dev(
		ldpc::ldpc_code_device(
			codeFile.c_str(), layerFile.c_str()));

	std::cout << "========================================================================================" << std::endl;
	code_dev->print();
	std::cout << "========================================================================================" << std::endl;

	if (code_dev->max_dc() > DEC_MAX_DC)
	{
		std::cout << "ERROR: maximum checknode degree(max_dc=" << code_dev->max_dc() << ") exceeds buffer(=" << DEC_MAX_DC << "). Adjust DEC_MAX_DC flag. Aborting...\n";
		exit(EXIT_FAILURE);
	}

	ldpc::cuda_ptr<ldpc::ldpc_sim_device> sim_dev(
		ldpc::ldpc_sim_device(
			code_dev, simFile.c_str(), mapFile.c_str(), numThreads));

	if (sim_dev->bits() > SIM_NUM_BITS)
	{
		std::cout << "ERROR: number of bits(bits=" << sim_dev->bits() << ") exceeds buffer(=" << SIM_NUM_BITS << "). Adjust SIM_NUM_BITS flag. Aborting...\n";
		exit(EXIT_FAILURE);
	}

	sim_dev->print();
	std::cout << "========================================================================================" << std::endl;

	//sim_dev->print();
#ifdef SHORT_LOG
	std::cout << "  FEC   |      FRAME     |   SNR   |    BER     |   FER      | AVGITERS  |  TIME/FRAME   \n";
	std::cout << "========+================+=========+============+============+===========+==============" << std::endl;
#endif

#if defined USE_LEGACY_DEC || defined USE_CPU_FRAME
	sim_dev->start();
#else
	sim_dev->start_device();
#endif

	return 0;
}
