#pragma once

#include <iostream>

namespace ldpc
{
    /**
     * @brief Implements Galois Field GF(2).
     * 
     */
    struct gf2
    {
        gf2() = default;
        gf2(const int val)
            : value(val != 0) {}

        gf2 &operator=(const int val)
        {
            value = (val != 0);
            return *this;
        }

        gf2 &operator+=(const gf2 &a)
        {
            value ^= a.value;
            return *this;
        }

        friend gf2 operator-(const gf2 &a);
        friend gf2 operator+(const gf2 &a, const gf2 &b);
        friend gf2 operator*(const gf2 &a, const gf2 &b);

        friend bool operator==(const gf2 &a, const gf2 &b);
        friend bool operator!=(const gf2 &a, const gf2 &b);

        friend std::ostream &operator<<(std::ostream &os, const gf2 &a);
        friend std::istream &operator>>(std::istream &is, gf2 &a);

        unsigned char value;
    };
} // namespace ldpc
