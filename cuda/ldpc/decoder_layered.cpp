#include "ldpc.h"
#include <exception>

using namespace ldpc;
using namespace std;


Ldpc_Decoder_cl::Ldpc_Decoder_cl(Ldpc_Code_cl* code) : ldpc_code(code)
{
    l_c2v = nullptr;
    l_v2c = nullptr;
    f = nullptr;
    b = nullptr;
    lsum = nullptr;

    c_out = nullptr;
    synd = nullptr;

    try
    {
        //num layers times num nnz
        l_c2v = new double[ldpc_code->nl() * ldpc_code->nnz()];
        l_v2c = new double[ldpc_code->nl() * ldpc_code->nnz()];
        f = new double[ldpc_code->nl() * ldpc_code->max_dc()];
        b = new double[ldpc_code->nl() * ldpc_code->max_dc()];

        lsum = new double[ldpc_code->nnz()];

        c_out = new bits_t[ldpc_code->nc()];
        synd = new bits_t[ldpc_code->nc()];
    }
    catch (exception& e)
    {
        cout << "Error: " << e.what() << endl;
        destroy_dec();
        exit(EXIT_FAILURE);
    }
}

Ldpc_Decoder_cl::~Ldpc_Decoder_cl() { destroy_dec(); }

void Ldpc_Decoder_cl::destroy_dec()
{
    if (l_c2v != nullptr)
        delete[] l_c2v;
    if (l_v2c != nullptr)
        delete[] l_v2c;
    if (f != nullptr)
        delete[] f;
    if (b != nullptr)
        delete[] b;
    if (lsum != nullptr)
        delete[] lsum;
    if (c_out != nullptr)
        delete[] c_out;
    if (synd != nullptr)
        delete[] synd;
}


const uint64_t Ldpc_Decoder_cl::decode_layered(double* llr_in, double* llr_out, const uint64_t& MaxIter, const uint8_t& early_termination)
{
    uint64_t I = 0;
    while (I < MaxIter)
    {
        /*
        //parallel
        for (size_t i = 0; i < code->nl; ++i)
            layered_dec(code, llr_in, l_c2v[i], lsum, l_v2c[i], code->layers[i], code->lw[i], f[i], b[i]);

        //interchange check node messages
        for(size_t i = 0; i < code->nnz; ++i)
        {
            lsum[i] = 0.0;
            for (size_t j = 0; j < code->nl; ++j)
                lsum[i] += l_c2v[j][i];
        }

        // app calculation
        for(size_t i = 0; i < code->nc; ++i)
        {
            llr_out[i] = llr_in[i];
            vn = code->vn[i];
            vw = code->vw[i];
            while(vw--)
                llr_out[i] += lsum[*vn++];
            c_out[i] = (llr_out[i] <= 0);
        }
        */

        ++I;

        if (early_termination)
        {
            if (isCodeword(c_out))
                break;
        }
    }

    return I;
}


const bool Ldpc_Decoder_cl::isCodeword(bits_t* c) //TODO - CUDA DEVICE FCT
{
    bool is_codeword = true;

    //calc syndrome
    bits_t s;
    for(size_t i = 0; i < ldpc_code->mc(); i++)
    {
        s = 0;
        for(size_t j = 0; j < ldpc_code->cw()[i]; j++)
            s ^= c[ldpc_code->c()[ldpc_code->cn()[i][j]]];

        if (s == 1)
        {
            is_codeword = false;
            break;
        }
    }

    return is_codeword;
}


/*
void layered_dec(ldpc_code_t* code, double* llr_in, double* l_c2v, double* l_c2v_sum, double* l_v2c, uint64_t* cn_subset, const uint64_t cn_size, double* f, double* b)
{
    size_t* vn;
    size_t* cn;

    size_t vw;
    size_t cw;


    for(size_t i = 0; i < code->nc; i++)
    {
        double tmp = llr_in[i];
        vw = code->vw[i];
        vn = code->vn[i];
        while(vw--)
            tmp += l_c2v_sum[*vn++];

        vn = code->vn[i];
        vw = code->vw[i];
        while(vw--)
        {
            l_v2c[*vn] = tmp - l_c2v[*vn];
            ++vn;
        }
    }


    for(size_t i = 0; i < cn_size; i++)
    {
        cw = code->cw[cn_subset[i]];
        cn = code->cn[cn_subset[i]];
        f[0] = l_v2c[*cn];
        b[cw-1] = l_v2c[*(cn+cw-1)];
        for(size_t j = 1; j < cw; j++)
        {
            f[j] = jacobian(f[j-1], l_v2c[*(cn+j)]);
            b[cw-1-j] = jacobian(b[cw-j], l_v2c[*(cn + cw-j-1)]);
        }

        l_c2v[*cn] = b[1];
        l_c2v[*(cn+cw-1)] = f[cw-2];

        for(size_t j = 1; j < cw-1; j++)
            l_c2v[*(cn+j)] = jacobian(f[j-1], b[j+1]);
    }
}
 */
