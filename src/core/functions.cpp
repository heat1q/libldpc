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
} // namespace ldpc
