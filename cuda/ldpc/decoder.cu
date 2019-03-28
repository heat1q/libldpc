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
	l_c2v_pre = nullptr;

	c_out = nullptr;

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

		l_c2v_pre = new double[num_layers * ldpc_code->nnz()]();
		lsum = new double[ldpc_code->nnz()]();

		c_out = new bits_t[ldpc_code->nc()]();
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
	l_c2v_pre = nullptr;

	c_out = nullptr;

	const uint64_t num_layers = ldpc_code->nl();

	//num layers times num nnz
	l_c2v = new double[num_layers * ldpc_code->nnz()]();
	l_v2c = new double[num_layers * ldpc_code->nnz()]();
	f = new double[num_layers * ldpc_code->max_dc()]();
	b = new double[num_layers * ldpc_code->max_dc()]();

	l_c2v_pre = new double[num_layers * ldpc_code->nnz()]();
	lsum = new double[ldpc_code->nnz()]();

	c_out = new bits_t[ldpc_code->nc()]();
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
	if (l_c2v_pre != nullptr)
		delete[] l_c2v_pre;
}


__host__ __device__ bool Ldpc_Decoder_cl::is_codeword()
{
    bool is_codeword = true;

    //calc syndrome
    bits_t s;
    for (size_t i = 0; i < ldpc_code->mc(); i++)
    {
        s = 0;
        for (size_t j = 0; j < ldpc_code->cw()[i]; j++)
            s ^= c_out[ldpc_code->c()[ldpc_code->cn()[i][j]]];

        if (s)
        {
            return false;
        }
    }

    return is_codeword;
}

uint64_t Ldpc_Decoder_cl::decode_legacy(double* llr_in, double* llr_out, const uint64_t& max_iter, const bool& early_termination)
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

        if (early_termination) {
            if (is_codeword()) {
                break;
            }
        }
    }

    return it;
}

uint64_t Ldpc_Decoder_cl::decode_layered_legacy(double* llr_in, double* llr_out, const uint64_t& max_iter, const bool& early_termination)
{
	size_t* vn;
    size_t* cn;

    size_t vw;
    size_t cw;

	size_t i_nnz;
	size_t i_dc;

    //initialize
    for (size_t i = 0; i < ldpc_code->nnz(); ++i)
    {
        lsum[i] = 0.0;
        for (size_t l = 0; l < ldpc_code->nl(); ++l)
        {
            l_c2v[l*ldpc_code->nnz()+i] = 0.0;
            l_v2c[l*ldpc_code->nnz()+i] = 0.0;
            l_c2v_pre[l*ldpc_code->nnz()+i] = 0.0;
        }
    }

    size_t I = 0;
    while (I < max_iter)
    {
        for (size_t l = 0; l < ldpc_code->nl(); ++l)
        {
			i_nnz = l*ldpc_code->nnz();
			i_dc = l*ldpc_code->max_dc();

            /* VN node intialization */
            for(size_t i = 0; i < ldpc_code->nc(); i++)
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
                    l_v2c[i_nnz + *vn] = tmp - l_c2v[i_nnz + *vn];
                    ++vn;
                }
            }

            //CN processing
            for(size_t i = 0; i < ldpc_code->lw()[l]; i++)
            {
                cw = ldpc_code->cw()[ldpc_code->layers()[l][i]];
                cn = ldpc_code->cn()[ldpc_code->layers()[l][i]];
                f[i_dc] = l_v2c[i_nnz + *cn];
                b[i_dc + cw-1] = l_v2c[i_nnz + *(cn+cw-1)];
                for(size_t j = 1; j < cw; j++)
                {
                    f[i_dc + j] = jacobian(f[i_dc + j-1], l_v2c[i_nnz + *(cn+j)]);
                    b[i_dc + cw-1-j] = jacobian(b[i_dc + cw-j], l_v2c[i_nnz + *(cn + cw-j-1)]);
                }

                l_c2v[i_nnz + *cn] = b[i_dc + 1];
                l_c2v[i_nnz + *(cn+cw-1)] = f[i_dc + cw-2];

                for(size_t j = 1; j < cw-1; j++)
                    l_c2v[i_nnz + *(cn+j)] = jacobian(f[i_dc + j-1], b[i_dc + j+1]);
            }

            //update the llr sum of layers, by replacing old llr of lyr l with new value
            for (size_t i = 0; i < ldpc_code->nnz(); ++i)
            {
                lsum[i] += l_c2v[i_nnz + i] - l_c2v_pre[i_nnz + i];
                l_c2v_pre[i_nnz + i] = l_c2v[i_nnz + i];
            }

            // app calculation
            for(size_t i = 0; i < ldpc_code->nc(); ++i)
            {
                llr_out[i] = llr_in[i];
                vn = ldpc_code->vn()[i];
                vw = ldpc_code->vw()[i];
                while(vw--)
                    llr_out[i] += lsum[*vn++];
                c_out[i] = (llr_out[i] <= 0);
            }

            if (early_termination)
            {
                if (is_codeword())
                    return I;
            }
        }

        ++I;
    }

    return I;
}


__global__ void cudakernel::setup_decoder(Ldpc_Code_cl* code_managed, Ldpc_Decoder_cl** dec_ptr)
{
	*dec_ptr = new Ldpc_Decoder_cl();
	(**dec_ptr).setup_decoder_device(code_managed);
	printf("Cuda Device :: Decoder set up!\n");
}


__global__ void cudakernel::destroy_decoder(Ldpc_Decoder_cl** dec_ptr)
{
	delete *dec_ptr;
	printf("Cuda Device :: Decoder destroyed!\n");
}


__global__ void cudakernel::clean_decoder(Ldpc_Decoder_cl** dec_ptr)
{
	uint_fast32_t index = blockIdx.x * blockDim.x + threadIdx.x;
	uint_fast32_t stride = blockDim.x * gridDim.x;

	for (size_t i = index; i < (**dec_ptr).ldpc_code->nnz(); i += stride)
	{
		(**dec_ptr).lsum[i] = 0.0;
		for (size_t l = 0; l < (**dec_ptr).ldpc_code->nl(); ++l)
		{
			(**dec_ptr).l_c2v[l*(**dec_ptr).ldpc_code->nnz()+i] = 0.0;
			(**dec_ptr).l_v2c[l*(**dec_ptr).ldpc_code->nnz()+i] = 0.0;
			(**dec_ptr).l_c2v_pre[l*(**dec_ptr).ldpc_code->nnz()+i] = 0.0;
		}
	}
}


__global__ void cudakernel::decode_lyr_vnupdate(Ldpc_Decoder_cl** dec_ptr, double* llr_in, size_t i_nnz)
{
	size_t* vn;
	size_t vw;

	uint_fast32_t index = blockIdx.x * blockDim.x + threadIdx.x;
	uint_fast32_t stride = blockDim.x * gridDim.x;

	//VN processing
	for (size_t i = index; i < (**dec_ptr).ldpc_code->nc(); i += stride)
	{
		double tmp = llr_in[i];
		vw = (**dec_ptr).ldpc_code->vw()[i];
		vn = (**dec_ptr).ldpc_code->vn()[i];
		while(vw--)
			tmp += (**dec_ptr).lsum[*vn++];

		vn = (**dec_ptr).ldpc_code->vn()[i];
		vw = (**dec_ptr).ldpc_code->vw()[i];
		while(vw--)
		{
			(**dec_ptr).l_v2c[i_nnz + *vn] = tmp - (**dec_ptr).l_c2v[i_nnz + *vn];
			++vn;
		}
	}
}


__global__ void cudakernel::decode_lyr_cnupdate(Ldpc_Decoder_cl** dec_ptr, size_t i_nnz, uint64_t L)
{
	size_t* cn;
	size_t cw;

	uint_fast32_t index = blockIdx.x * blockDim.x + threadIdx.x;
	uint_fast32_t stride = blockDim.x * gridDim.x;

	double f_tmp[10];//[sizeof((**dec_ptr).f[0])/8] = {0}; // - TODO
	double b_tmp[10];//[sizeof((**dec_ptr).b[0])/8] = {0}; // - TODO

	//CN processing
	for (size_t i = index; i < (**dec_ptr).ldpc_code->lw()[L]; i += stride)
	{
		cw = (**dec_ptr).ldpc_code->cw()[(**dec_ptr).ldpc_code->layers()[L][i]];
		cn = (**dec_ptr).ldpc_code->cn()[(**dec_ptr).ldpc_code->layers()[L][i]];
		f_tmp[0] = (**dec_ptr).l_v2c[i_nnz + *cn];
		b_tmp[cw-1] = (**dec_ptr).l_v2c[i_nnz + *(cn+cw-1)];
		for(size_t j = 1; j < cw; j++)
		{
			f_tmp[j] = jacobian(f_tmp[j-1], (**dec_ptr).l_v2c[i_nnz + *(cn+j)]);
			b_tmp[cw-1-j] = jacobian(b_tmp[cw-j], (**dec_ptr).l_v2c[i_nnz + *(cn + cw-j-1)]);
		}

		(**dec_ptr).l_c2v[i_nnz + *cn] = b_tmp[1];
		(**dec_ptr).l_c2v[i_nnz + *(cn+cw-1)] = f_tmp[cw-2];

		for(size_t j = 1; j < cw-1; j++)
			(**dec_ptr).l_c2v[i_nnz + *(cn+j)] = jacobian(f_tmp[j-1], b_tmp[j+1]);
	}
}


__global__ void cudakernel::decode_lyr_sumllr(Ldpc_Decoder_cl** dec_ptr, size_t i_nnz)
{
	uint_fast32_t index = blockIdx.x * blockDim.x + threadIdx.x;
	uint_fast32_t stride = blockDim.x * gridDim.x;

	//sum llrs
	for (size_t i = index; i < (**dec_ptr).ldpc_code->nnz(); i += stride)
	{
		(**dec_ptr).lsum[i] += (**dec_ptr).l_c2v[i_nnz + i] - (**dec_ptr).l_c2v_pre[i_nnz + i];
		(**dec_ptr).l_c2v_pre[i_nnz + i] = (**dec_ptr).l_c2v[i_nnz + i];
	}
}


__global__ void cudakernel::decode_lyr_appcalc(Ldpc_Decoder_cl** dec_ptr, double* llr_in, double* llr_out)
{
	size_t* vn;
	size_t vw;

	uint_fast32_t index = blockIdx.x * blockDim.x + threadIdx.x;
	uint_fast32_t stride = blockDim.x * gridDim.x;

	//app calc
	for (size_t i = index; i < (**dec_ptr).ldpc_code->nc(); i += stride)
	{
		llr_out[i] = llr_in[i];
		vn = (**dec_ptr).ldpc_code->vn()[i];
		vw = (**dec_ptr).ldpc_code->vw()[i];
		while(vw--)
			llr_out[i] += (**dec_ptr).lsum[*vn++];
		(**dec_ptr).c_out[i] = (llr_out[i] <= 0);
	}
}


//tmpl fcts need definition in each file?
template<typename T> void ldpc::printVector(T *x, const size_t &l)
{
    cout << "[";
    for (size_t i = 0; i < l-1; ++i)
        cout << x[i] << " ";
    cout << x[l-1] << "]";
}
