#include "ldpc.cuh"
#include <exception>

using namespace ldpc;
using namespace std;


Ldpc_Decoder_cl::Ldpc_Decoder_cl(Ldpc_Code_cl* code, const uint16_t I, const bool early_term, const bool mgd)
	: ldpc_code(code), max_iter(I), early_termination(early_term)
{
	if (mgd)
	{
		setup_dec_mgd();
		prefetch_dec();
	}
	else
	{
		//setup_dec();
	}

	isMgd = mgd;
}


Ldpc_Decoder_cl::~Ldpc_Decoder_cl()
{
	if (isMgd)
	{
		destroy_dec_mgd();
	}
	else
	{
		//destroy_dec();
	}
}


void Ldpc_Decoder_cl::setup_dec_mgd()
{
	block_size = 256;
    num_blocks = ceil((ldpc_code->nnz() + block_size - 1) / block_size);

	l_c2v = nullptr;
	l_v2c = nullptr;
	f = nullptr;
	b = nullptr;
	fb_ref = nullptr;
	lsum = nullptr;
	l_c2v_pre = nullptr;
	c_out = nullptr;
	synd = nullptr;
	llr_in = nullptr;
	llr_out = nullptr;

	const uint64_t num_layers = ldpc_code->nl();

	try
	{
		//num layers times num nnz
		cudaMallocManaged(&l_c2v, sizeof(double)*num_layers*ldpc_code->nnz());
		if (l_c2v == NULL || l_c2v == nullptr) { throw runtime_error("l_c2v alloc failed."); }
		cudaMallocManaged(&l_v2c, sizeof(double)*num_layers*ldpc_code->nnz());
		if (l_v2c == NULL || l_v2c == nullptr) { throw runtime_error("l_v2c alloc failed."); }
		cudaMallocManaged(&l_c2v_pre, sizeof(double)*num_layers*ldpc_code->nnz());
		if (l_c2v_pre == NULL || l_c2v_pre == nullptr) { throw runtime_error("l_c2v_pre alloc failed."); }

		cudaMallocManaged(&f, sizeof(double)*num_layers*ldpc_code->max_dc());
		if (f == NULL || f == nullptr) { throw runtime_error("f alloc failed."); }
		cudaMallocManaged(&b, sizeof(double)*num_layers*ldpc_code->max_dc());
		if (b == NULL || b == nullptr) { throw runtime_error("b alloc failed."); }

		cudaMallocManaged(&fb_ref, ldpc_code->max_dc());
		if (fb_ref == NULL || fb_ref == nullptr) { throw runtime_error("fb_ref alloc failed."); }

		cudaMallocManaged(&lsum, sizeof(double)*ldpc_code->nnz());
		if (lsum == NULL || lsum == nullptr) { throw runtime_error("lsum alloc failed."); }

		cudaMallocManaged(&llr_in, sizeof(double)*ldpc_code->nc());
		if (llr_in == NULL || llr_in == nullptr) { throw runtime_error("llr_in alloc failed."); }
		cudaMallocManaged(&llr_out, sizeof(double)*ldpc_code->nc());
		if (llr_out == NULL || llr_out == nullptr) { throw runtime_error("llr_out alloc failed."); }
		cudaMallocManaged(&c_out, sizeof(bits_t)*ldpc_code->nc());
		if (c_out == NULL || c_out == nullptr) { throw runtime_error("c_out alloc failed."); }

		cudaMallocManaged(&synd, sizeof(bits_t)*ldpc_code->mc());
		if (synd == NULL || synd == nullptr) { throw runtime_error("synd alloc failed."); }
	}
	catch (exception& e)
	{
		cout << "Error: " << e.what() << endl;
		destroy_dec_mgd();
		exit(EXIT_FAILURE);
	}
}

void Ldpc_Decoder_cl::prefetch_dec()
{
	cudaDeviceSynchronize();

	int dev = -1;
	cudaGetDevice(&dev);

	const uint64_t num_layers = ldpc_code->nl();

	cudaMemPrefetchAsync(l_c2v, sizeof(double)*num_layers*ldpc_code->nnz(), dev, NULL);
	cudaMemPrefetchAsync(l_v2c, sizeof(double)*num_layers*ldpc_code->nnz(), dev, NULL);
	cudaMemPrefetchAsync(l_c2v_pre, sizeof(double)*num_layers*ldpc_code->nnz(), dev, NULL);

	cudaMemPrefetchAsync(fb_ref, ldpc_code->max_dc(), dev, NULL);

	cudaMemPrefetchAsync(lsum, sizeof(double)*ldpc_code->nnz(), dev, NULL);

	cudaMemPrefetchAsync(llr_in, sizeof(double)*ldpc_code->nc(), dev, NULL);
	cudaMemPrefetchAsync(llr_out, sizeof(double)*ldpc_code->nc(), dev, NULL);
	cudaMemPrefetchAsync(c_out, sizeof(double)*ldpc_code->nc(), dev, NULL);

	cudaMemPrefetchAsync(synd, sizeof(double)*ldpc_code->mc(), dev, NULL);

	cudaMemPrefetchAsync(this, sizeof(Ldpc_Decoder_cl), dev, NULL);
}

void Ldpc_Decoder_cl::destroy_dec_mgd()
{
	if (l_c2v != nullptr) { cudaFree(l_c2v); }
	if (l_v2c != nullptr) { cudaFree(l_v2c); }
	if (l_c2v_pre != nullptr) { cudaFree(l_c2v_pre); }
	if (f != nullptr) { cudaFree(f); }
	if (b != nullptr) { cudaFree(b); }
	if (fb_ref != nullptr) { cudaFree(fb_ref); }
	if (lsum != nullptr) { cudaFree(lsum); }
	if (c_out != nullptr) { cudaFree(c_out); }
	if (synd != nullptr) { cudaFree(synd); }
	if (llr_in != nullptr) { cudaFree(llr_in); }
	if (llr_out != nullptr) { cudaFree(llr_out); }
}


__host__ __device__ bool Ldpc_Decoder_cl::is_codeword()
{
	is_cw = true;

	//calc syndrome
	cudakernel::decoder::calc_synd<<<get_num_size(ldpc_code->mc(), 256), 256>>>(this);
	cudaDeviceSynchronize();

	return is_cw;
}


__host__ __device__ bool Ldpc_Decoder_cl::is_codeword_legacy()
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


uint16_t Ldpc_Decoder_cl::decode_legacy()
{
    uint16_t it;

    size_t* vn;
    size_t* cn;

    size_t vw;
    size_t cw;

    //initialize with llrs
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

        // VN node processing
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
            if (is_codeword_legacy()) {
                break;
            }
        }
    }

    return it;
}


uint16_t Ldpc_Decoder_cl::decode_layered_legacy()
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

    uint16_t I = 0;
    while (I < max_iter)
    {
        for (size_t l = 0; l < ldpc_code->nl(); ++l)
        {
			i_nnz = l*ldpc_code->nnz();
			i_dc = l*ldpc_code->max_dc();

            // VN node intialization
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
                if (is_codeword_legacy())
				{
                    return I;
				}
            }
        }

        ++I;
    }

    return I;
}


uint16_t Ldpc_Decoder_cl::decode_layered()
{
	cudakernel::decoder::decode_layered<<<1, 1>>>(this);
	cudaDeviceSynchronize();

	return iter;
}


//tmpl fcts need definition in each file?
template<typename T> void ldpc::printVector(T *x, const size_t &l)
{
    cout << "[";
    for (size_t i = 0; i < l-1; ++i)
        cout << x[i] << " ";
    cout << x[l-1] << "]";
}


/*
	Cudakernels
*/
__global__ void cudakernel::decoder::clean_decoder(Ldpc_Decoder_cl* dec_mgd)
{
	uint_fast32_t index = blockIdx.x * blockDim.x + threadIdx.x;
	uint_fast32_t stride = blockDim.x * gridDim.x;

	for (size_t i = index; i < dec_mgd->ldpc_code->nnz(); i += stride)
	{
		dec_mgd->lsum[i] = 0.0;
		for (size_t l = 0; l < dec_mgd->ldpc_code->nl(); ++l)
		{
			dec_mgd->l_c2v[l*dec_mgd->ldpc_code->nnz()+i] = 0.0;
			dec_mgd->l_v2c[l*dec_mgd->ldpc_code->nnz()+i] = 0.0;
			dec_mgd->l_c2v_pre[l*dec_mgd->ldpc_code->nnz()+i] = 0.0;
		}
	}
}


__global__ void cudakernel::decoder::decode_layered(Ldpc_Decoder_cl* dec_mgd)
{
    size_t i_nnz;

    //zero everything out
    cudakernel::decoder::clean_decoder<<<get_num_size(dec_mgd->ldpc_code->nnz(), 256), 256>>>(dec_mgd);

    uint16_t I = 0;
    for (; I < dec_mgd->max_iter; ++I)
    {
        for (uint64_t l = 0; l < dec_mgd->ldpc_code->nl(); ++l)
        {
            i_nnz = dec_mgd->ldpc_code->nnz()*l;

            //launch kernels here
            cudakernel::decoder::decode_lyr_vnupdate<<<get_num_size(dec_mgd->ldpc_code->nc(), 256), 256>>>(dec_mgd, i_nnz);
            cudakernel::decoder::decode_lyr_cnupdate<<<get_num_size(dec_mgd->ldpc_code->lw()[l], 256), 256>>>(dec_mgd, i_nnz, l);
            cudakernel::decoder::decode_lyr_sumllr<<<get_num_size(dec_mgd->ldpc_code->nnz(), 256), 256>>>(dec_mgd, i_nnz);
            cudakernel::decoder::decode_lyr_appcalc<<<get_num_size(dec_mgd->ldpc_code->nc(), 256), 256>>>(dec_mgd);

            if (dec_mgd->early_termination)
            {
                if (dec_mgd->is_codeword_legacy()) //break
                {
                    l = dec_mgd->ldpc_code->nl();
					I += dec_mgd->max_iter;
                }
            }
        }
    }

    cudaDeviceSynchronize();
	dec_mgd->iter = I % dec_mgd->max_iter;
}


__global__ void cudakernel::decoder::decode_lyr_vnupdate(Ldpc_Decoder_cl* dec_mgd, size_t i_nnz)
{
	size_t* vn;
	size_t vw;

	uint_fast32_t index = blockIdx.x * blockDim.x + threadIdx.x;
	uint_fast32_t stride = blockDim.x * gridDim.x;

	//VN processing
	for (size_t i = index; i < dec_mgd->ldpc_code->nc(); i += stride)
	{
		double tmp = dec_mgd->llr_in[i];
		vw =  dec_mgd->ldpc_code->vw()[i];
		vn = dec_mgd->ldpc_code->vn()[i];
		while(vw--)
			tmp += dec_mgd->lsum[*vn++];

		vn = dec_mgd->ldpc_code->vn()[i];
		vw = dec_mgd->ldpc_code->vw()[i];
		while(vw--)
		{
			dec_mgd->l_v2c[i_nnz + *vn] = tmp - dec_mgd->l_c2v[i_nnz + *vn];
			++vn;
		}
	}
}


__global__ void cudakernel::decoder::decode_lyr_cnupdate(Ldpc_Decoder_cl* dec_mgd, size_t i_nnz, uint64_t l)
{
	size_t* cn;
	size_t cw;

	uint_fast32_t index = blockIdx.x * blockDim.x + threadIdx.x;
	uint_fast32_t stride = blockDim.x * gridDim.x;

	double f_tmp[sizeof(dec_mgd->fb_ref)];
	double b_tmp[sizeof(dec_mgd->fb_ref)];

	//CN processing
	for (size_t i = index; i < dec_mgd->ldpc_code->lw()[l]; i += stride)
	{
		cw = dec_mgd->ldpc_code->cw()[dec_mgd->ldpc_code->layers()[l][i]];
		cn = dec_mgd->ldpc_code->cn()[dec_mgd->ldpc_code->layers()[l][i]];
		f_tmp[0] = dec_mgd->l_v2c[i_nnz + *cn];
		b_tmp[cw-1] = dec_mgd->l_v2c[i_nnz + *(cn+cw-1)];
		for(size_t j = 1; j < cw; j++)
		{
			f_tmp[j] = jacobian(f_tmp[j-1], dec_mgd->l_v2c[i_nnz + *(cn+j)]);
			b_tmp[cw-1-j] = jacobian(b_tmp[cw-j], dec_mgd->l_v2c[i_nnz + *(cn + cw-j-1)]);
		}

		dec_mgd->l_c2v[i_nnz + *cn] = b_tmp[1];
		dec_mgd->l_c2v[i_nnz + *(cn+cw-1)] = f_tmp[cw-2];

		for(size_t j = 1; j < cw-1; j++)
			dec_mgd->l_c2v[i_nnz + *(cn+j)] = jacobian(f_tmp[j-1], b_tmp[j+1]);
	}
}


__global__ void cudakernel::decoder::decode_lyr_sumllr(Ldpc_Decoder_cl* dec_mgd, size_t i_nnz)
{
	uint_fast32_t index = blockIdx.x * blockDim.x + threadIdx.x;
	uint_fast32_t stride = blockDim.x * gridDim.x;

	//sum llrs
	for (size_t i = index; i < dec_mgd->ldpc_code->nnz(); i += stride)
	{
		dec_mgd->lsum[i] += dec_mgd->l_c2v[i_nnz + i] - dec_mgd->l_c2v_pre[i_nnz + i];
		dec_mgd->l_c2v_pre[i_nnz + i] = dec_mgd->l_c2v[i_nnz + i];
	}
}


__global__ void cudakernel::decoder::decode_lyr_appcalc(Ldpc_Decoder_cl* dec_mgd)
{
	size_t* vn;
	size_t vw;

	uint_fast32_t index = blockIdx.x * blockDim.x + threadIdx.x;
	uint_fast32_t stride = blockDim.x * gridDim.x;

	//app calc
	for (size_t i = index; i < dec_mgd->ldpc_code->nc(); i += stride)
	{
		dec_mgd->llr_out[i] = dec_mgd->llr_in[i];
		vn = dec_mgd->ldpc_code->vn()[i];
		vw = dec_mgd->ldpc_code->vw()[i];
		while(vw--)
			dec_mgd->llr_out[i] += dec_mgd->lsum[*vn++];
		dec_mgd->c_out[i] = (dec_mgd->llr_out[i] <= 0);
	}
}


__global__ void cudakernel::decoder::calc_synd(Ldpc_Decoder_cl* dec_mgd)
{
	uint_fast32_t index = blockIdx.x * blockDim.x + threadIdx.x;
	uint_fast32_t stride = blockDim.x * gridDim.x;

	for (size_t i = index; i < dec_mgd->ldpc_code->mc(); i += stride)
	{
		dec_mgd->synd[i] = 0;
		for (size_t j = 0; j < dec_mgd->ldpc_code->cw()[i]; j++)
		{
			dec_mgd->synd[i] ^= dec_mgd->c_out[dec_mgd->ldpc_code->c()[dec_mgd->ldpc_code->cn()[i][j]]];
		}

		if (dec_mgd->synd[i])
		{
			dec_mgd->is_cw = false;
		}
	}
}
