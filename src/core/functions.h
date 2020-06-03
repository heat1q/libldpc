#pragma once

#include <iostream>
#include <cinttypes>
#include <vector>
#include <memory>
#include <cmath>
#include <cstring>
#include <chrono>
#include <fstream>
#include <random>

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

namespace pgd
{
using bits_t = std::uint8_t;
using labels_t = std::uint16_t;
using symbols_t = std::uint16_t;

using vec_bits_t = std::vector<bits_t>;
using vec_labels_t = std::vector<labels_t>;
using vec_symbols_t = std::vector<symbols_t>;
using vec_size_t = std::vector<std::size_t>;
using vec_double_t = std::vector<double>;

using mat_bits_t = std::vector<std::vector<bits_t>>;
using mat_size_t = std::vector<std::vector<std::size_t>>;
using mat_double_t = std::vector<std::vector<double>>;

void dec2bin(std::size_t val, uint8_t m);
int sign(double a);

} // namespace pgd