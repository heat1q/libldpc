#include "ldpcsim.h"

using namespace ldpc;

// /usr/local/cuda-9.2/bin/nvcc -x cu -std=c++11 sim_cuda.cpp ldpc/ldpc.cpp ldpc/decoder.cpp device/kernel.cpp device/ldpcsimdevice.cpp  -o sim_cuda -arch sm_35 -rdc=true -O3 -w
int main(int argc, char *argv[])
{
	bool abort = false;

	std::string codeFile;
	std::string simFile;
	std::string mapFile;
	std::string layerFile;
	uint16_t numThreads;

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
				numThreads = static_cast<unsigned short>(atoi(argv[2 * i + 1]));
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

	if (abort)
	{
		std::cout << "======================== LDPC Simulation ========================\n";
		std::cout << "                        (Layered Decoding)                       \n";
		std::cout << "                         Usage Reminder:                         \n";
		std::cout << "         Main -code CodeFile -sim SimFile -map MappingFile       \n";
		std::cout << "               optional: -layer LayerFile                        \n";
		std::cout << "                         -threads NumThreads                     \n";
		std::cout << "                                                                 \n";
		std::cout << "                                                                 \n";
		std::cout << "                CodeFile: Name of the code file                  \n";
		std::cout << "               CodeMapFile: Name of mapping file                 \n";
		std::cout << "               SimFile: Name of simulation file                  \n";
		std::cout << "              LayerFile: Name of code layer file                 \n";
		std::cout << "                         for layered decoding                    \n";
		std::cout << "              NumThreads: Number of threads for parallel         \n";
		std::cout << "                          processing of frames (default: 1)      \n";
		std::cout << "=================================================================" << std::endl;
		exit(EXIT_FAILURE);
	}

	cuda_ptr<ldpc_code_device> code_dev(
		ldpc_code_device(
			codeFile.c_str(), layerFile.c_str()));

	cuda_ptr<ldpc_sim_device> sim_dev(
		ldpc_sim_device(
			code_dev, simFile.c_str(), mapFile.c_str()));

	sim_dev->print();
	//sim_dev->start();
	sim_dev->start_device();

	return 0;
}
