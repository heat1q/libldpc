#include "functions.h"

namespace pgd
{
/**
 * @brief Prints a decimal as binary with m bits
 * 
 * @param val 
 * @param m 
 */
void dec2bin(std::size_t val, uint8_t m)
{
    for (std::size_t i = 0; i < m; i++)
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
} // namespace pgd