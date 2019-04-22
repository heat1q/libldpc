#include "../ldpcsim.h"

using namespace ldpc;

/*
__global__ void cudakernel::sim::awgn()
{
	printf("curand\n");

	curandState_t state[10];
	for (size_t i = 0; i < 10; i++) {
		curand_init(clock64(), 1, 0, &state[i]);
	}

	double a[10];
	for (size_t i = 0; i < 10; i++) {
		a[i] = curand_normal(&state[i]);
		printf("%.3f\n", a[i]);
	}
}
*/
__global__ void cudakernel::sim::frame_proc(cudamgd_ptr<ldpc_sim_device> pSim)
{
	
}

__global__ void cudakernel::sim::awgn(cudamgd_ptr<ldpc_sim_device> pSim, double sigma2)
{
	const size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t sx = blockDim.x * gridDim.x;

	double a = 0;

	for (size_t i = ix; i < pSim->n(); i += sx)
	{
		a = curand_normal(&(pSim->mCurandState[i])) * sqrt(sigma2);
		pSim->mY[i] = pSim->cstll().X()[pSim->mX[i]] + a;
	}
}

__global__ void cudakernel::sim::setup_randn(cudamgd_ptr<ldpc_sim_device> pSim)
{
	// initialize curand
	const size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t sx = blockDim.x * gridDim.x;

	for (size_t i = ix; i < pSim->n(); i += sx)
	{
		curand_init(clock64(), i, 0, &(pSim->mCurandState[i]));
	}
}

/*
 *	Decoder kernels
 */
__global__ void cudakernel::decoder::clean_decoder(ldpc_decoder_device* pDecMgd)
{
	const size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t sx = blockDim.x * gridDim.x;

    for (size_t i = ix; i < pDecMgd->mLSum.size(); i += sx)
    {
        pDecMgd->mLSum[i] = 0.0;
        for (size_t l = 0; l < pDecMgd->mLdpcCode->nl(); ++l)
        {
            pDecMgd->mLc2v[l*pDecMgd->mLdpcCode->nnz()+i] = 0.0;
            pDecMgd->mLv2c[l*pDecMgd->mLdpcCode->nnz()+i] = 0.0;
            pDecMgd->mLc2vPre[l*pDecMgd->mLdpcCode->nnz()+i] = 0.0;
        }
    }
}


__global__ void cudakernel::decoder::decode_layered(ldpc_decoder_device* pDecMgd)
{
    size_t pI;

	const size_t gridSizeNC = get_num_size(pDecMgd->mLdpcCode->nc(), NUM_THREADS);
	const size_t gridSizeNNZ = get_num_size(pDecMgd->mLdpcCode->nnz(), NUM_THREADS);

    //zero everything out
    cudakernel::decoder::clean_decoder<<<gridSizeNNZ, NUM_THREADS>>>(pDecMgd);

    uint16_t I = 0;
    while (I < pDecMgd->max_iter())
    {
        for (uint64_t l = 0; l < pDecMgd->mLdpcCode->nl(); ++l)
        {
            pI = pDecMgd->mLdpcCode->nnz()*l;

            //launching kernels
            cudakernel::decoder::decode_lyr_vnupdate<<<gridSizeNC, NUM_THREADS>>>(pDecMgd, pI);
            cudakernel::decoder::decode_lyr_cnupdate<<<get_num_size(pDecMgd->mLdpcCode->lw()[l], NUM_THREADS/2), NUM_THREADS/2>>>(pDecMgd, pI, l);
            cudakernel::decoder::decode_lyr_sumllr<<<gridSizeNNZ, NUM_THREADS>>>(pDecMgd, pI);
            cudakernel::decoder::decode_lyr_appcalc<<<gridSizeNC, NUM_THREADS>>>(pDecMgd);

            if (pDecMgd->early_termination())
            {
                if (pDecMgd->is_codeword()) //break
                {
                    goto break_here;
                }
            }
        }

        ++I;
    }

    break_here:
    cudaDeviceSynchronize();

    pDecMgd->mIter = I;
}


__global__ void cudakernel::decoder::decode_lyr_vnupdate(ldpc_decoder_device* pDecMgd, size_t pI)
{
    size_t* vn;
    size_t vw;

    const size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t sx = blockDim.x * gridDim.x;

    //VN processing
    for (size_t i = ix; i < pDecMgd->mLdpcCode->nc(); i += sx)
    {
        double tmp = pDecMgd->mLLRIn[i];
        vw = pDecMgd->mLdpcCode->vn()[i].size();
        vn = pDecMgd->mLdpcCode->vn()[i].data();
        while(vw--)
            tmp += pDecMgd->mLSum[*vn++];

        vn = pDecMgd->mLdpcCode->vn()[i].data();
        vw = pDecMgd->mLdpcCode->vn()[i].size();
        while(vw--)
        {
            pDecMgd->mLv2c[pI + *vn] = tmp - pDecMgd->mLc2v[pI + *vn];
            ++vn;
        }
    }
}


__global__ void cudakernel::decoder::decode_lyr_cnupdate(ldpc_decoder_device* pDecMgd, size_t pI, uint64_t pL)
{
    size_t* cn;
    size_t cw;

    const size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t sx = blockDim.x * gridDim.x;

    double f_tmp[30];
    double b_tmp[30];

    //CN processing
    for (size_t i = ix; i < pDecMgd->mLdpcCode->layers()[pL].size(); i += sx)
    {
        cw = pDecMgd->mLdpcCode->cn()[pDecMgd->mLdpcCode->layers()[pL][i]].size();
        cn = pDecMgd->mLdpcCode->cn()[pDecMgd->mLdpcCode->layers()[pL][i]].data();
        f_tmp[0] = pDecMgd->mLv2c[pI + *cn];
        b_tmp[cw-1] = pDecMgd->mLv2c[pI + *(cn+cw-1)];
        for(size_t j = 1; j < cw; j++)
        {
            f_tmp[j] = jacobian(f_tmp[j-1], pDecMgd->mLv2c[pI + *(cn+j)]);
            b_tmp[cw-1-j] = jacobian(b_tmp[cw-j], pDecMgd->mLv2c[pI + *(cn + cw-j-1)]);
        }

        pDecMgd->mLc2v[pI + *cn] = b_tmp[1];
        pDecMgd->mLc2v[pI + *(cn+cw-1)] = f_tmp[cw-2];

        for(size_t j = 1; j < cw-1; j++) {
            pDecMgd->mLc2v[pI + *(cn+j)] = jacobian(f_tmp[j-1], b_tmp[j+1]);
		}
    }
}


__global__ void cudakernel::decoder::decode_lyr_sumllr(ldpc_decoder_device* pDecMgd, size_t pI)
{
	const size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t sx = blockDim.x * gridDim.x;

    //sum llrs
    for (size_t i = ix; i < pDecMgd->mLdpcCode->nnz(); i += sx)
    {
        pDecMgd->mLSum[i] += pDecMgd->mLc2v[pI + i] - pDecMgd->mLc2vPre[pI + i];
        pDecMgd->mLc2vPre[pI + i] = pDecMgd->mLc2v[pI + i];
    }
}


__global__ void cudakernel::decoder::decode_lyr_appcalc(ldpc_decoder_device* pDecMgd)
{
    size_t* vn;
    size_t vw;

	const size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t sx = blockDim.x * gridDim.x;

    //app calc
    for (size_t i = ix; i < pDecMgd->mLdpcCode->nc(); i += sx)
    {
        pDecMgd->mLLROut[i] = pDecMgd->mLLRIn[i];
        vn = pDecMgd->mLdpcCode->vn()[i].data();
        vw = pDecMgd->mLdpcCode->vn()[i].size();
        while(vw--)
            pDecMgd->mLLROut[i] += pDecMgd->mLSum[*vn++];
        pDecMgd->mCO[i] = (pDecMgd->mLLROut[i] <= 0);
    }
}


__global__ void cudakernel::decoder::calc_synd(ldpc_decoder_device* pDecMgd)
{
	const size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t sx = blockDim.x * gridDim.x;

    for (size_t i = ix; i < pDecMgd->mLdpcCode->mc(); i += sx)
    {
        pDecMgd->mSynd[i] = 0;
        for (auto cni : pDecMgd->mLdpcCode->cn()[i])
        {
			pDecMgd->mSynd[i] ^= pDecMgd->mCO[pDecMgd->mLdpcCode->c()[cni]];
        }

        if (pDecMgd->mSynd[i])
        {
            pDecMgd->mIsCW = false;
        }
    }
}

/*
__global__ void cudakernel::decoder::clean_decoder(ldpc_decoder* pDecMgd)
{
    const size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t sx = blockDim.x * gridDim.x;

    for (size_t i = ix; i < pDecMgd->mLdpcCode->nnz(); i += sx)
    {
        pDecMgd->mLSum[i] = 0.0;
        for (size_t l = 0; l < pDecMgd->mLdpcCode->nl(); ++l)
        {
            pDecMgd->mLc2v[l*pDecMgd->mLdpcCode->nnz()+i] = 0.0;
            pDecMgd->mLv2c[l*pDecMgd->mLdpcCode->nnz()+i] = 0.0;
            pDecMgd->mLc2vPre[l*pDecMgd->mLdpcCode->nnz()+i] = 0.0;
        }
    }
}


__global__ void cudakernel::decoder::decode_layered(ldpc_decoder* pDecMgd)
{
    size_t pI;

    const size_t gridSizeNC = get_num_size(pDecMgd->mLdpcCode->nc(), NUM_THREADS);
    const size_t gridSizeNNZ = get_num_size(pDecMgd->mLdpcCode->nnz(), NUM_THREADS);

    //zero everything out
    cudakernel::decoder::clean_decoder<<<gridSizeNNZ, NUM_THREADS>>>(pDecMgd);

    uint16_t I = 0;
    while (I < pDecMgd->max_iter())
    {
        for (uint64_t l = 0; l < pDecMgd->mLdpcCode->nl(); ++l)
        {
            pI = pDecMgd->mLdpcCode->nnz()*l;

            //launching kernels
            cudakernel::decoder::decode_lyr_vnupdate<<<gridSizeNC, NUM_THREADS>>>(pDecMgd, pI);
            cudakernel::decoder::decode_lyr_cnupdate<<<get_num_size(pDecMgd->mLdpcCode->lw()[l], NUM_THREADS/2), NUM_THREADS/2>>>(pDecMgd, pI, l);
            cudakernel::decoder::decode_lyr_sumllr<<<gridSizeNNZ, NUM_THREADS>>>(pDecMgd, pI);
            cudakernel::decoder::decode_lyr_appcalc<<<gridSizeNC, NUM_THREADS>>>(pDecMgd);

            if (pDecMgd->early_termination())
            {
                if (pDecMgd->is_codeword()) //break
                {
                    //l = pDecMgd->mLdpcCode->nl();
                    //I += pDecMgd->max_iter();
                    goto break_here;
                }
            }
        }

        ++I;
    }

break_here:
    cudaDeviceSynchronize();

    pDecMgd->mIter = I;
}


__global__ void cudakernel::decoder::decode_lyr_vnupdate(ldpc_decoder* pDecMgd, size_t pI)
{
    size_t* vn;
    size_t vw;

    const size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t sx = blockDim.x * gridDim.x;

    //VN processing
    for (size_t i = ix; i < pDecMgd->mLdpcCode->nc(); i += sx)
    {
        double tmp = pDecMgd->mLLRIn[i];
        vw =  pDecMgd->mLdpcCode->vw()[i];
        vn = pDecMgd->mLdpcCode->vn()[i];
        while(vw--)
            tmp += pDecMgd->mLSum[*vn++];

        vn = pDecMgd->mLdpcCode->vn()[i];
        vw = pDecMgd->mLdpcCode->vw()[i];
        while(vw--)
        {
            pDecMgd->mLv2c[pI + *vn] = tmp - pDecMgd->mLc2v[pI + *vn];
            ++vn;
        }
    }
}


__global__ void cudakernel::decoder::decode_lyr_cnupdate(ldpc_decoder* pDecMgd, size_t pI, uint64_t pL)
{
    size_t* cn;
    size_t cw;

    const size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t sx = blockDim.x * gridDim.x;

    double f_tmp[sizeof(pDecMgd->FBREF)];
    double b_tmp[sizeof(pDecMgd->FBREF)];

    //CN processing
    for (size_t i = ix; i < pDecMgd->mLdpcCode->lw()[pL]; i += sx)
    {
        cw = pDecMgd->mLdpcCode->cw()[pDecMgd->mLdpcCode->layers()[pL][i]];
        cn = pDecMgd->mLdpcCode->cn()[pDecMgd->mLdpcCode->layers()[pL][i]];
        f_tmp[0] = pDecMgd->mLv2c[pI + *cn];
        b_tmp[cw-1] = pDecMgd->mLv2c[pI + *(cn+cw-1)];
        for(size_t j = 1; j < cw; j++)
        {
            f_tmp[j] = jacobian(f_tmp[j-1], pDecMgd->mLv2c[pI + *(cn+j)]);
            b_tmp[cw-1-j] = jacobian(b_tmp[cw-j], pDecMgd->mLv2c[pI + *(cn + cw-j-1)]);
        }

        pDecMgd->mLc2v[pI + *cn] = b_tmp[1];
        pDecMgd->mLc2v[pI + *(cn+cw-1)] = f_tmp[cw-2];

        for(size_t j = 1; j < cw-1; j++) {
            pDecMgd->mLc2v[pI + *(cn+j)] = jacobian(f_tmp[j-1], b_tmp[j+1]);
        }
    }
}


__global__ void cudakernel::decoder::decode_lyr_sumllr(ldpc_decoder* pDecMgd, size_t pI)
{
    const size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t sx = blockDim.x * gridDim.x;

    //sum llrs
    for (size_t i = ix; i < pDecMgd->mLdpcCode->nnz(); i += sx)
    {
        pDecMgd->mLSum[i] += pDecMgd->mLc2v[pI + i] - pDecMgd->mLc2vPre[pI + i];
        pDecMgd->mLc2vPre[pI + i] = pDecMgd->mLc2v[pI + i];
    }
}


__global__ void cudakernel::decoder::decode_lyr_appcalc(ldpc_decoder* pDecMgd)
{
    size_t* vn;
    size_t vw;

    const size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t sx = blockDim.x * gridDim.x;

    //app calc
    for (size_t i = ix; i < pDecMgd->mLdpcCode->nc(); i += sx)
    {
        pDecMgd->mLLROut[i] = pDecMgd->mLLRIn[i];
        vn = pDecMgd->mLdpcCode->vn()[i];
        vw = pDecMgd->mLdpcCode->vw()[i];
        while(vw--)
            pDecMgd->mLLROut[i] += pDecMgd->mLSum[*vn++];
        pDecMgd->mCO[i] = (pDecMgd->mLLROut[i] <= 0);
    }
}


__global__ void cudakernel::decoder::calc_synd(ldpc_decoder* pDecMgd)
{
    const size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t sx = blockDim.x * gridDim.x;

    for (size_t i = ix; i < pDecMgd->mLdpcCode->mc(); i += sx)
    {
        pDecMgd->mSynd[i] = 0;
        for (size_t j = 0; j < pDecMgd->mLdpcCode->cw()[i]; j++)
        {
            pDecMgd->mSynd[i] ^= pDecMgd->mCO[pDecMgd->mLdpcCode->c()[pDecMgd->mLdpcCode->cn()[i][j]]];
        }

        if (pDecMgd->mSynd[i])
        {
            pDecMgd->mIsCW = false;
        }
    }
}

__global__ void cudakernel::sim::sim_calc_llrs(ldpc_sim_device* pSim, double pSigma2)
{
	const size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t sx = blockDim.x * gridDim.x;

}


*/

/*
__global__ void cudakernel::sim::sim_test(Ldpc_Decoder_cl* dec_mgd)
{
	curandState_t state;
	curand_init(clock64(), 1, 0, &state);
	for (size_t i=0; i<dec_mgd->ldpc_code->nc(); ++i)
	{
		dec_mgd->llr_in[i] = curand_normal(&state);
		dec_mgd->llr_out[i] = 0.0;
	}
}
*/
