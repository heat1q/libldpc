#include "../src/core/ldpc.h"

namespace ldpc_tests
{
    void gf2()
    {
        try
        {
            auto check = ldpc::bits_t(0);
            if ((check * 1) != 0) throw 0;
            if ((check + 1) == 0) throw 0;
            if ((check + 1 + 1) != 0) throw 0;
            if ((-check) == 0) throw 0;
        }
        catch (...)
        {
            throw std::runtime_error("failed: gf2 arithmetics");
        }
        std::cout << "passed: gf2 arithmetics" << std::endl;
    }

    void is_generator_matrix(const ldpc::ldpc_code &code)
    {
        for (int i = 0; i < code.mc(); i++) // rows of H
        {
            for (int j = 0; j < code.kc(); j++) // cols of G transpose = rows of G
            {
                auto &h = code.H().row_neighbor()[i];
                auto &g = code.G().row_neighbor()[j];

                auto sum = ldpc::bits_t(0);
                for (auto hi : h)
                {
                    for (auto gj : g)
                    {
                        sum += (hi.nodeIndex == gj.nodeIndex) *
                               code.G().nz_entry()[gj.edgeIndex].value *
                               code.H().nz_entry()[hi.edgeIndex].value;
                    }
                }

                if (sum != 0)
                {
                    throw std::runtime_error("failed: is_generator_matrix");
                }
            }
        }

        std::cout << "passed: is_generator_matrix" << std::endl;
    }

    void codeword(const ldpc::ldpc_code &code)
    {
        ldpc::vec_bits_t u(code.kc(), 1);
        for (auto &x : u)
        {
            x = rand() % 2;
        }
        auto cw = code.G().multiply_left(u);
        auto synd = code.H().multiply_right(cw);

        for (auto s : synd)
        {
            if (s != 0)
            {
                throw std::runtime_error("failed: encoding random information word");
            }
        }

        std::cout << "passed: encoding random information word" << std::endl;
    }
} // namespace ldpc_tests
