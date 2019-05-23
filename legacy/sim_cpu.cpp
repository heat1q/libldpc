#include "sim/ldpcsim.h"


int main(int argc, char *argv[])
{
	bool abort = false;

	std::string codeFile;
	std::string simFile;
	std::string mapFile;
	std::string layerFile;

if (argc == 7 || argc == 9)
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


	//some blabla
	std::cout << "==================================== LDPC Simulation ===================================\n";
	std::cout << "codefile: " << codeFile << "\n";
	std::cout << "simfile: " << simFile << "\n";
	std::cout << "mappingfile: " << mapFile << "\n";
	if (!layerFile.empty()) 
	{ 
		std::cout << "layerfile: " << layerFile << "\n";
	}
	
	std::cout << "\nDFLAGS: \t";
#ifdef LOG_FRAME_TIME
	std::cout << "LOG_FRAME_TIME ";
#endif
#ifdef LOG_TP
	std::cout << "LOG_TP ";
#endif
#ifdef CN_APPROX_LIN
	std::cout << "CN_APPROX_LIN ";
#endif
#ifdef CN_APPROX_MINSUM
	std::cout << "CN_APPROX_MINSUM ";
#endif

	std::cout << "\n";
	if (abort)
	{
		std::cout << "==================================== USAGE REMINDER ====================================\n";
		std::cout << "                    ./Main  -code CodeFile -sim SimFile -map MappingFile                 \n";
		std::cout << "                   Optinal: -layer LayerFile                                             \n";
		std::cout << "                                                                                         \n";
		std::cout << "                  CodeFile: Name of the code file                                        \n";
		std::cout << "               CodeMapFile: Name of mapping file                                         \n";
		std::cout << "                   SimFile: Name of simulation file                                      \n";
		std::cout << "                 LayerFile: Name of layer file                                           \n";
		std::cout << "========================================================================================" << std::endl;
		exit(EXIT_FAILURE);
	}

	ldpc::ldpc_code code(codeFile.c_str(), layerFile.c_str());
	ldpc::ldpc_decoder decoder(&code, 0, true);

	std::cout << "========================================================================================" << std::endl;
	code.print();
	std::cout << "========================================================================================" << std::endl;

	ldpc::ldpc_sim sim(&code, &decoder, simFile.c_str(), mapFile.c_str());

	sim.print();

#ifdef LOG_TP
	std::cout << "===================================================================================================================" << std::endl;
	std::cout << "  FEC   |      FRAME     |   SNR   |    BER     |    FER     | AVGITERS  |   T FRAME   |   T DEC    |  THROUGHPUT  \n";
	std::cout << "========+================+=========+============+============+===========+=============+============+==============" << std::endl;
#else
	std::cout << "========================================================================================" << std::endl;
	std::cout << "  FEC   |      FRAME     |   SNR   |    BER     |    FER     | AVGITERS  |  TIME/FRAME   \n";
	std::cout << "========+================+=========+============+============+===========+==============" << std::endl;
#endif

	sim.start();

	return 0;
}
