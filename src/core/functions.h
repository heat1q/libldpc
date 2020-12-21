#pragma once

#include <cinttypes>
#include <vector>
#include <memory>
#include <cmath>
#include <cstring>
#include <chrono>
#include <fstream>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <variant>
#include <forward_list>

#include "gf2.h"
#include "sparse.h"

#define TIME_PROF(log, exec, unit)                                                                                                             \
    do                                                                                                                                         \
    {                                                                                                                                          \
        std::string str_unit = std::string(unit);                                                                                              \
        float a = 1;                                                                                                                           \
        if (str_unit == std::string("s"))                                                                                                      \
        {                                                                                                                                      \
            a = 1e-9;                                                                                                                          \
        }                                                                                                                                      \
        else if (str_unit == std::string("ms"))                                                                                                \
        {                                                                                                                                      \
            a = 1e-6;                                                                                                                          \
        }                                                                                                                                      \
        else if (str_unit == std::string("us"))                                                                                                \
        {                                                                                                                                      \
            a = 1e-3;                                                                                                                          \
        }                                                                                                                                      \
        else if (str_unit == std::string("ns"))                                                                                                \
        {                                                                                                                                      \
            a = 1;                                                                                                                             \
        }                                                                                                                                      \
        else                                                                                                                                   \
        {                                                                                                                                      \
            a = 1;                                                                                                                             \
            str_unit = std::string("ns");                                                                                                      \
        }                                                                                                                                      \
        auto start = std::chrono::high_resolution_clock::now();                                                                                \
        exec;                                                                                                                                  \
        auto elapsed = std::chrono::high_resolution_clock::now() - start;                                                                      \
        printf("[TIMEPROF]: " log ": ");                                                                                                       \
        printf("%.3f %s\n", static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed).count()) * a, str_unit.c_str()); \
    } while (0);

namespace ldpc
{
    /**
    * @brief Define outstream operator for vector to ease debugging.
    * 
    * @tparam T 
    * @param os 
    * @param v 
    * @return std::ostream& 
    */
    template <class T>
    std::ostream &operator<<(std::ostream &os, const std::vector<T> &v)
    {
        os << "[";
        typename std::vector<T>::const_iterator it = v.begin();
        for (; it != v.end(); ++it)
        {
            os << *it;
            if (it != v.end() - 1)
            {
                os << ", ";
            }
        }
        os << "]";
        return os;
    }

    template <class T>
    std::ostream &operator<<(std::ostream &os, const std::forward_list<T> &l)
    {
        for (const auto &x : l)
        {
            os << x << " -> ";
        }
        os << "NULL";
        return os;
    }

    using bits_t = gf2;
    using u64 = unsigned long;
    using u32 = unsigned int;
    using u8 = unsigned char;

    using vec_bits_t = std::vector<bits_t>;
    using vec_u64 = std::vector<u64>;
    using vec_double_t = std::vector<double>;
    using vec_int = std::vector<int>;

    using mat_bits_t = std::vector<std::vector<bits_t>>;
    using mat_u64 = std::vector<std::vector<u64>>;
    using mat_double_t = std::vector<std::vector<double>>;
    using mat_int = std::vector<std::vector<int>>;

    constexpr u8 ERASURE = 'E';

    struct
    {
        bool earlyTerm;
        u32 iterations;
        const char *type;
    } typedef decoder_param;

    struct
    {
        u64 seed;
        double xRange[3];
        const char *type;
    } typedef channel_param;

    struct
    {
        u32 threads;
        u64 maxFrames;
        u64 fec;
        const char *resultFile;
    } typedef simulation_param;

    std::ostream &operator<<(std::ostream &os, const decoder_param &p);
    std::ostream &operator<<(std::ostream &os, const channel_param &p);
    std::ostream &operator<<(std::ostream &os, const simulation_param &p);

    void dec2bin(u64 val, uint8_t m);

} // namespace ldpc
