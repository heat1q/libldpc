#include "functions.h"

namespace ldpc
{
    /**
    * @brief Prints a decimal as binary with m bits
    * 
    * @param val 
    * @param m 
    */
    void dec2bin(u64 val, uint8_t m)
    {
        for (u64 i = 0; i < m; i++)
        {
            printf("%lu", (val >> (m - i - 1) & 0x01));
        }
    }

    /**
    * @brief Sign function
    * 
    * @param a 
    * @return int 
    */
    int sign(double a)
    {
        return (a <= 0) ? -1 : 1;
    }

    std::ostream &operator<<(std::ostream &os, const decoder_param &p)
    {
        os << "Type: " << p.type << "\n";
        os << "Iterations: " << p.iterations << "\n";
        os << "Early Termination: " << p.earlyTerm;
        return os;
    }

    std::ostream &operator<<(std::ostream &os, const channel_param &p)
    {
        os << "Type: " << p.type << "\n";
        os << "Seed: " << p.seed << "\n";
        os << "Range: " << p.xRange << "\n";
        os << "Values: " << p.xVals;
        return os;
    }

    std::ostream &operator<<(std::ostream &os, const simulation_param &p)
    {
        os << "Threads: " << p.threads << "\n";
        os << "FEC: " << p.fec << "\n";
        os << "Max Frames: " << p.maxFrames << "\n";
        os << "Output File: " << p.resultFile;
        return os;
    }
} // namespace ldpc
