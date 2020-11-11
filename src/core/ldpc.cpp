#include "ldpc.h"
#include <iterator>

namespace ldpc
{

    ldpc_code::ldpc_code(const std::string &pcFileName)
        : mMaxDC(0),
          mH(),
          mG()
    {
        try
        {
            read_H(pcFileName);
        }
        catch (std::exception &e)
        {
            std::cout << "Error: ldpc_code(): " << e.what() << std::endl;
            exit(EXIT_FAILURE);
        }
    }

    ldpc_code::ldpc_code(const std::string &pcFileName, const std::string &genFileName)
        : ldpc_code(pcFileName)
    {
        if (!genFileName.empty())
        {
            try
            {
                read_G(genFileName);
            }
            catch (std::exception &e)
            {
                std::cout << "Error: ldpc_code(): " << e.what() << std::endl;
                exit(EXIT_FAILURE);
            }
        }
    }

    void ldpc_code::read_H(const std::string &pcFileName)
    {
        std::ifstream infile(pcFileName);
        std::string line;
        int skipLines = 0;

        if (!infile.good()) throw std::runtime_error("can not open file for reading");

        while (getline(infile, line))
        {
            // support legacy code files
            // where code parameters are separated by ':'
            if (auto i = line.find(':'); i != std::string::npos)
            {
                int index;
                auto token = line.substr(0, i);
                std::istringstream record(line.substr(i+1));

                if (token.find("puncture") != std::string::npos)
                {
                    while (record >> index) mPuncture.push_back(index);
                } 
                else if (token.find("shorten") != std::string::npos)
                {
                    while (record >> index) mShorten.push_back(index);
                }

                ++skipLines;
            }
            else
            {
                break;
            }
        }

        infile.close();

        mH.read_from_file(pcFileName, skipLines);

        // maximum check node degree
        auto tmp = std::max_element(mH.row_neighbor().begin(), mH.row_neighbor().end(), [](const auto &a, const auto &b) { return (a.size() < b.size()); });
        mMaxDC = tmp->size();

        // position of transmitted bits
        for (int i = 0; i < nc(); i++)
        {
            auto tmp = std::find(mShorten.cbegin(), mShorten.cend(), i);
            if (tmp != mShorten.cend()) continue; // skip if current index shortened
            tmp = std::find(mPuncture.cbegin(), mPuncture.cend(), i);
            if (tmp != mPuncture.cend()) continue; // skip if current index punctured

            mBitPos.push_back(i);
        }
    }

    void ldpc_code::read_G(const std::string &genFileName)
    {
        mG.read_from_file(genFileName, 0);
    }

    /**
    * @brief Prints parameters of LDPC code
    * 
    */
    std::ostream &operator<<(std::ostream &os, const ldpc_code &code)
    {
        // calculate real rate of transmitted code
        auto rate = 1. - static_cast<double>(code.mct()) / static_cast<double>(code.nct());

        os << "N : " << code.nc() << "\n";
        os << "M : " << code.mc() << "\n";
        os << "K : " << code.kc() << "\n";
        os << "NNZ : " << code.nnz() << "\n";
        //os << "Rank: " << code.mRank << "\n";
        //os << "max dc : " << code.max_dc() << "\n";
        os << "puncture[" << code.puncture().size() << "] : " << code.puncture() << "\n";
        os << "shorten[" << code.shorten().size() << "] : " << code.shorten() << "\n";
        os << "Rate : " << rate << "\n";
        os << "N (transmitted) : " << code.nct() << "\n";
        os << "M (transmitted) : " << code.mct() << "\n";
        os << "K (transmitted) : " << code.kct() << "\n";
        return os;
    }
} // namespace ldpc
