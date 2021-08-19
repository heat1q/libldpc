#include "sim/ldpcsim.h"
#include "../include/argparse/argparse.hpp"


int main(int argc, char *argv[])
{
    argparse::ArgumentParser parser("ldpc");
    parser.add_argument("codefile").help("LDPC parity-check matrix file containing all non-zero entries.");
    parser.add_argument("output-file").help("Results output file.");
    parser.add_argument("snr-range").help("{MIN} {MAX} {STEP}").nargs(3).action([](const std::string &s) { return std::stod(s); });

    parser.add_argument("-G", "--gen-matrix").help("Generator matrix file.").default_value(std::string(""));

    parser.add_argument("-i", "--num-iterations").help("Number of iterations for decoding. (Default: 50)").default_value(unsigned(50)).action([](const std::string &s) { return static_cast<unsigned>(std::stoul(s)); });
    parser.add_argument("-s", "--seed").help("RNG seed. (Default: 0)").default_value(ldpc::u64(0)).action([](const std::string &s) { return std::stoul(s); });
    parser.add_argument("-t", "--num-threads").help("Number of frames to be decoded in parallel. (Default: 1)").default_value(unsigned(1)).action([](const std::string &s) { return static_cast<unsigned>(std::stoul(s)); });

    parser.add_argument("--channel").help("Specifies channel: \"AWGN\", \"BSC\", \"BEC\" (Default: AWGN)").default_value(std::string("AWGN"));
    parser.add_argument("--decoding").help("Specifies decoding algorithm: \"BP\", \"BP_MS\" (Default: BP)").default_value(std::string("BP"));
    parser.add_argument("--max-frames").help("Limit number of decoded frames.").default_value(ldpc::u64(10e9)).action([](const std::string &s) { return std::stoul(s); });
    parser.add_argument("--frame-error-count").help("Maximum frame errors for given simulation point.").default_value(ldpc::u64(50)).action([](const std::string &s) { return std::stoul(s); });
    parser.add_argument("--no-early-term").help("Disable early termination for decoding.").default_value(false).implicit_value(true);

    try
    {
        parser.parse_args(argc, argv);

        auto snr = parser.get<ldpc::vec_double_t>("snr-range");
        if (snr[0] > snr[1]) throw std::runtime_error("snr min > snr max");
        
        auto code = std::make_shared<ldpc::ldpc_code>(parser.get<std::string>("codefile"), parser.get<std::string>("-G"));
        std::cout << "========================================================================================" << std::endl;
        std::cout << "Parity-Check Matrix: " << parser.get<std::string>("codefile") << std::endl;
        std::cout << "Generator Matrix: " << parser.get<std::string>("-G") << std::endl;
        std::cout << *code << std::endl;
        std::cout << "========================================================================================" << std::endl;

        std::string decType, chType, resFile;
        decType = parser.get<std::string>("--decoding");
        chType = parser.get<std::string>("--channel");
        resFile = parser.get<std::string>("output-file");

        // decoder parameters
        ldpc::decoder_param decoderParams;
        decoderParams.iterations = parser.get<ldpc::u32>("--num-iterations");
        decoderParams.earlyTerm = !parser.get<bool>("--no-early-term");
        decoderParams.type = decType.c_str();

        // channel parameters
        ldpc::channel_param channelParams;
        channelParams.seed = parser.get<ldpc::u64>("--seed");
        channelParams.type = chType.c_str();
        for (int i = 0; i < 3; ++i)
        {
           channelParams.xRange[i] = snr[i];
        }

        // simulation parameters
        ldpc::simulation_param simulationParams;
        simulationParams.threads = parser.get<ldpc::u32>("--num-threads");
        simulationParams.fec = parser.get<ldpc::u64>("--frame-error-count");
        simulationParams.maxFrames = parser.get<ldpc::u64>("--max-frames");
        simulationParams.resultFile = resFile.c_str();

        ldpc::ldpc_sim sim(
            code,
            decoderParams,
            channelParams,
            simulationParams
        );
        std::cout << sim << std::endl;
        std::cout << "========================================================================================" << std::endl;

        bool stop = false;
        sim.start(&stop);
    }
    catch (const std::runtime_error &e)
    {
        std::cout << e.what() << std::endl;
        std::cout << parser;
        exit(EXIT_FAILURE);
    }

    return 0;
}
