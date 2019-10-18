#include "sim/ldpcsim.h"

int main(int argc, char *argv[])
{
    bool abort = false;

    std::string codeFile;
    std::string simFile;
    std::string mapFile("");
    std::uint16_t numThreads = 1;

    if (argc == 5 || argc == 7 || argc == 9)
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
            else if (strcmp(argv[2 * i + 1], "-threads") == 0)
            {
                numThreads = std::stoul(argv[2 * i + 2]);
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
    std::cout << "threads: " << numThreads << "\n";

    std::cout << "\nDFLAGS: \t";
#ifdef LOG_FRAME_TIME
    std::cout << "LOG_FRAME_TIME ";
#endif
#ifdef CN_APPROX_MINSUM
    std::cout << "CN_APPROX_MINSUM ";
#endif
#ifdef LOG_CW
    std::cout << "LOG_CW ";
#endif

    std::cout << "\n";
    if (abort)
    {
        std::cout << "==================================== USAGE REMINDER ====================================\n";
        std::cout << "                    ./Main  -code CodeFile -sim SimFile -map MappingFile                 \n";
        std::cout << "                  optional: -map MappingFile=''                                          \n";
        std::cout << "                            -threads NumThreads=1                                        \n";
        std::cout << "                                                                                         \n";
        std::cout << "                  CodeFile: Name of the code file                                        \n";
        std::cout << "                   SimFile: Name of simulation file                                      \n";
        std::cout << "               MappingFile: Name of mapping file                                         \n";
        std::cout << "                NumThreads: Number of parallel threads                                   \n";
        std::cout << "========================================================================================" << std::endl;
        exit(EXIT_FAILURE);
    }

    ldpc::ldpc_code code(codeFile.c_str());

    std::cout << "========================================================================================" << std::endl;
    code.print();
    std::cout << "========================================================================================" << std::endl;

    ldpc::ldpc_sim sim(&code, simFile.c_str(), mapFile.c_str(), numThreads);

    sim.print();

    std::cout << "========================================================================================" << std::endl;
    std::cout << "  FEC   |      FRAME     |   SNR   |    BER     |    FER     | AVGITERS  |  TIME/FRAME   \n";
    std::cout << "========+================+=========+============+============+===========+==============" << std::endl;

    sim.start();

    return 0;
}
