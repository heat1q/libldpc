#include "gf2.h"

namespace ldpc
{
    gf2 operator-(const gf2 &a)
    {
        return gf2(~a.value);
    }

    gf2 operator+(const gf2 &a, const gf2 &b)
    {
        return gf2(a.value ^ b.value);
    }

    gf2 operator*(const gf2 &a, const gf2 &b)
    {
        return gf2(a.value & b.value);
    }

    bool operator==(const gf2 &a, const gf2 &b)
    {
        return (a.value == b.value);
    }

    bool operator!=(const gf2 &a, const gf2 &b)
    {
        return (a.value != b.value);
    }

    std::ostream &operator<<(std::ostream &os, const gf2 &a)
    {
        os << static_cast<int>(a.value);
        return os;
    }

    std::istream &operator>>(std::istream &is, gf2 &a)
    {
        int val;
        is >> val;
        a.value = val != 0;
        return is;
    }
} // namespace ldpc
