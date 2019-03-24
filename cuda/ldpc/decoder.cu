#include "ldpc.cuh"
#include <exception>

using namespace ldpc;
using namespace std;

__host__ __device__ Ldpc_Decoder_cl::Ldpc_Decoder_cl() {}
Ldpc_Decoder_cl::Ldpc_Decoder_cl(Ldpc_Code_cl* code) { setup_decoder(code); }
__host__ __device__ Ldpc_Decoder_cl::~Ldpc_Decoder_cl()
{
	if (init)
		destroy_dec();
}

void Ldpc_Decoder_cl::setup_decoder(Ldpc_Code_cl* code)
{
	init = true;
	ldpc_code = code;

	l_c2v = nullptr;
	l_v2c = nullptr;
	f = nullptr;
	b = nullptr;
	lsum = nullptr;

	c_out = nullptr;
	synd = nullptr;

	#ifdef QC_LYR_DEC
	const uint64_t num_layers = ldpc_code->nl();
	#else
	const uint64_t num_layers = 1;
	#endif

	try
	{
		//num layers times num nnz
		l_c2v = new double[num_layers * ldpc_code->nnz()]();
		l_v2c = new double[num_layers * ldpc_code->nnz()]();
		f = new double[num_layers * ldpc_code->max_dc()]();
		b = new double[num_layers * ldpc_code->max_dc()]();

		lsum = new double[ldpc_code->nnz()]();

		c_out = new bits_t[ldpc_code->nc()]();
		synd = new bits_t[ldpc_code->nc()]();
	}
	catch (exception& e)
	{
		cout << "Error: " << e.what() << endl;
		destroy_dec();
		exit(EXIT_FAILURE);
	}
}

__device__ void Ldpc_Decoder_cl::setup_decoder_device(Ldpc_Code_cl* code)
{
	init = true;
	ldpc_code = code;

	l_c2v = nullptr;
	l_v2c = nullptr;
	f = nullptr;
	b = nullptr;
	lsum = nullptr;

	c_out = nullptr;
	synd = nullptr;

	const uint64_t num_layers = ldpc_code->nl();

	//num layers times num nnz
	l_c2v = new double[num_layers * ldpc_code->nnz()]();
	l_v2c = new double[num_layers * ldpc_code->nnz()]();
	f = new double[num_layers * ldpc_code->max_dc()]();
	b = new double[num_layers * ldpc_code->max_dc()]();

	lsum = new double[ldpc_code->nnz()]();

	c_out = new bits_t[ldpc_code->nc()]();
	synd = new bits_t[ldpc_code->nc()]();
}

__host__ __device__ void Ldpc_Decoder_cl::destroy_dec()
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

#ifdef QC_LYR_DEC
uint64_t Ldpc_Decoder_cl::decode_layered(double* llr_in, double* llr_out, const uint64_t& MaxIter, const bool& early_termination)
{
    uint64_t I = 0;
    while (I < MaxIter)
    {
        decode_lyr_nodeupdate_global(llr_in);
        decode_lyr_sumllr_global();
        decode_lyr_appcalc_global(llr_in, llr_out);

        ++I;

        if (early_termination)
        {
            if (is_codeword_global(c_out))
                break;
        }
    }

    return I;
}

void Ldpc_Decoder_cl::decode_lyr_nodeupdate_global(double* llr_in) //TODO - CUDA DEVICE FCT
{
    size_t index_msg;
    size_t index_fb;

    size_t* cn_subset;

    size_t* vn;
    size_t* cn;

    size_t vw;
    size_t cw;

    for (size_t L=0; L<ldpc_code->nl(); ++L) //parallelize
    {
        index_msg = L*ldpc_code->nnz();
        index_fb = L*ldpc_code->max_dc();

        cn_subset = ldpc_code->layers()[L];

        //VN init
        for (size_t i = 0; i < ldpc_code->nc(); i++)
        {
            double tmp = llr_in[i];
            vw = ldpc_code->vw()[i];
            vn = ldpc_code->vn()[i];
            while(vw--)
                tmp += lsum[*vn++];

            vn = ldpc_code->vn()[i];
            vw = ldpc_code->vw()[i];
            while(vw--)
            {
                l_v2c[index_msg + *vn] = tmp - l_c2v[index_msg + *vn];
                ++vn;
            }
        }

        //CN update
        for (size_t i = 0; i < ldpc_code->lw()[L]; i++)
        {
            cw = ldpc_code->cw()[cn_subset[i]];
            cn = ldpc_code->cn()[cn_subset[i]];
            f[index_fb] = l_v2c[index_msg + *cn];
            b[index_fb + cw-1] = l_v2c[index_msg + *(cn+cw-1)];
            for(size_t j = 1; j < cw; j++)
            {
                f[index_fb + j] = jacobian(f[index_fb + j-1], l_v2c[index_msg + *(cn+j)]);
                b[index_fb + cw-1-j] = jacobian(b[index_fb + cw-j], l_v2c[index_msg + *(cn + cw-j-1)]);
            }

            l_c2v[index_msg + *cn] = b[index_fb + 1];
            l_c2v[index_msg + *(cn+cw-1)] = f[index_fb + cw-2];

            for(size_t j = 1; j < cw-1; j++)
                l_c2v[index_msg + *(cn+j)] = jacobian(f[index_fb + j-1], b[index_fb + j+1]);
        }
    }
}

void Ldpc_Decoder_cl::decode_lyr_sumllr_global() //TODO - CUDA DEVICE FCT
{
    for(size_t i = 0; i < ldpc_code->nnz(); ++i)
    {
        lsum[i] = 0.0;
        for (size_t L = 0; L < ldpc_code->nl(); ++L)
            lsum[i] += l_c2v[L*ldpc_code->nnz() + i];
    }
}

void Ldpc_Decoder_cl::decode_lyr_appcalc_global(double* llr_in, double* llr_out) //TODO - CUDA DEVICE FCT
{
    size_t* vn;
    size_t vw;
    for(size_t i = 0; i < ldpc_code->nc(); ++i)
    {
        llr_out[i] = llr_in[i];
        vn = ldpc_code->vn()[i];
        vw = ldpc_code->vw()[i];
        while(vw--)
            llr_out[i] += lsum[*vn++];
        c_out[i] = (llr_out[i] <= 0);
    }
}
#endif

bool Ldpc_Decoder_cl::is_codeword_global(bits_t* c) //TODO - CUDA DEVICE FCT
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

uint64_t Ldpc_Decoder_cl::decode(double* llr_in, double* llr_out, const uint64_t& max_iter, const bool& early_termination)
{
    size_t it;

    size_t* vn;
    size_t* cn;

    size_t vw;
    size_t cw;

    /* initialize with llrs */
    for(size_t i = 0; i < ldpc_code->nnz(); i++) {
        l_v2c[i] = llr_in[ldpc_code->c()[i]];
    }

    it = 0;
    while(it < max_iter) {
        for(size_t i = 0; i < ldpc_code->mc(); i++) {
            cw = ldpc_code->cw()[i];
            cn = ldpc_code->cn()[i];
            f[0] = l_v2c[*cn];
            b[cw-1] = l_v2c[*(cn+cw-1)];
            for(size_t j = 1; j < cw; j++) {
                f[j] = jacobian(f[j-1], l_v2c[*(cn+j)]);
                b[cw-1-j] = jacobian(b[cw-j], l_v2c[*(cn + cw-j-1)]);
            }

            l_c2v[*cn] = b[1];
            l_c2v[*(cn+cw-1)] = f[cw-2];
            for(size_t j = 1; j < cw-1; j++) {
                l_c2v[*(cn+j)] = jacobian(f[j-1], b[j+1]);
            }
        }

        /* VN node processing */
        for(size_t i = 0; i < ldpc_code->nc(); i++) {
            double tmp = llr_in[i];
            vw = ldpc_code->vw()[i];
            vn = ldpc_code->vn()[i];
            while(vw--) {
                tmp += l_c2v[*vn++];
            }
            vn = ldpc_code->vn()[i];
            vw = ldpc_code->vw()[i];
            while(vw--) {
                l_v2c[*vn] = tmp - l_c2v[*vn];
                vn++;
            }
        }

        // app calculation
        for(size_t i = 0; i < ldpc_code->nc(); i++) {
            llr_out[i] = llr_in[i];
            vn = ldpc_code->vn()[i];
            vw = ldpc_code->vw()[i];
            while(vw--) {
                llr_out[i] += l_c2v[*vn++];
            }
            c_out[i] = (llr_out[i] <= 0);
        }

        it++;

        if(early_termination) {
            if(is_codeword_global(c_out)) {
                break;
            }
        }
    }

    return it;
}

//tmpl fcts need definition in each file?
template<typename T> void ldpc::printVector(T *x, const size_t &l)
{
    cout << "[";
    for (size_t i = 0; i < l-1; ++i)
        cout << x[i] << " ";
    cout << x[l-1] << "]";
}
