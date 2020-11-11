#include "../src/core/ldpc.h"
#include "../include/argparse/argparse.hpp"

#include "ldpctest.cpp"

int main(int argc, char *argv[])
{
    argparse::ArgumentParser parser("ldpc_tests");
    parser.add_argument("codefile").help("LDPC codefile containing all non-zero entries, compressed sparse row (CSR) format.");
    parser.add_argument("-G").help("Generator matrix,  compressed sparse row (CSR) format.").default_value(std::string(""));

    std::string pcFile, genFile;
    try
    {
        parser.parse_args(argc, argv);
        pcFile = parser.get<std::string>("codefile");
        genFile = parser.get<std::string>("-G");
    }
    catch (const std::runtime_error &e)
    {
        std::cout << e.what() << std::endl;
        std::cout << parser;
        exit(EXIT_FAILURE);
    }

    try
    {
        ldpc::ldpc_code code(pcFile, genFile);
        ldpc_tests::gf2();
        ldpc_tests::rank(code);
        ldpc_tests::is_generator_matrix(code);
        ldpc_tests::codeword(code);

        std::cout << "All tests passed." << std::endl;
    }
    catch (const std::exception &e)
    {
        std::cout << "Assessment failed: " << e.what() << std::endl;
    }

    return 0;
}
