#include "../ldpcsim.h"

using namespace ldpc;

/*
*	Simulation kernels
*/
__global__ void cudakernel::sim::setup_rng(ldpc_sim_device *pSim)
{
    // initialize curand
    const size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t sx = blockDim.x * gridDim.x;

    for (size_t i = ix; i < pSim->n(); i += sx)
    {
        curand_init(clock64(), i, 0, &(pSim->mCurandState[i]));
    }

    for (size_t j = ix; j < pSim->ldpc_code()->nct(); j += sx)
    {
        curand_init(clock64(), j, 0, &(pSim->mCurandStateEncoding[j]));
    }
}

__global__ void cudakernel::sim::frame_proc(ldpc_sim_device *pSim, double pSigma2)
{
    //encodeall0
    cudakernel::sim::encode_all0<<<get_num_size(pSim->ldpc_code()->nct(), NUM_THREADS), NUM_THREADS>>>(pSim);

    //map c to x
    cudakernel::sim::map_c_to_x<<<get_num_size(pSim->n(), NUM_THREADS), NUM_THREADS>>>(pSim);

    //sim awgn
    cudakernel::sim::awgn<<<get_num_size(pSim->n(), NUM_THREADS), NUM_THREADS>>>(pSim, pSigma2);

    //num & punture
    //TODO

    //calc llrin
    cudakernel::sim::calc_llrs<<<get_num_size(pSim->n(), NUM_THREADS), NUM_THREADS>>>(pSim, pSigma2);
    cudakernel::sim::calc_llrin<<<get_num_size(pSim->ldpc_code()->nc(), NUM_THREADS), NUM_THREADS>>>(pSim);

    //decode
    cudakernel::decoder::decode_layered<<<1, 1>>>(pSim->mLdpcDecoder.get());

    cudaDeviceSynchronize();
}

//encode the input to all zero
__global__ void cudakernel::sim::encode_all0(ldpc_sim_device *pSim)
{
    const size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t sx = blockDim.x * gridDim.x;

    for (size_t i = ix; i < pSim->ldpc_code()->nct(); i += sx)
    {
        pSim->mC[pSim->bits_pos()[i]] = (curand_uniform(&(pSim->mCurandStateEncoding[i])) > 0.5000);
    }

    for (size_t i = ix; i < pSim->ldpc_code()->num_puncture(); i += sx)
    {
        pSim->mC[pSim->ldpc_code()->puncture()[i]] = (curand_uniform(&(pSim->mCurandStateEncoding[i])) > 0.5000);
    }

    for (size_t i = 0; i < pSim->ldpc_code()->num_shorten(); i++)
    {
        pSim->mC[pSim->ldpc_code()->shorten()[i]] = 0;
    }

    //map_c_to_x()
}

__global__ void cudakernel::sim::awgn(ldpc_sim_device *pSim, double pSigma2)
{
    const size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t sx = blockDim.x * gridDim.x;

    double a = 0;

    for (size_t i = ix; i < pSim->n(); i += sx)
    {
        a = curand_normal(&(pSim->mCurandState[i])) * sqrt(pSigma2);
        pSim->mY[i] = pSim->cstll().X()[pSim->mX[i]] + a;
        //Pn += a * a;
        //Px += pSim->cstll().X()[pSim->mX[i]] * pSim->cstll().X()[pSim->mX[i]];
    }

    //return Px/Pn
}

__global__ void cudakernel::sim::calc_llrs(ldpc_sim_device *pSim, double pSigma2)
{
    const size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t sx = blockDim.x * gridDim.x;

    double llr_tmp[SIM_NUM_BITS];

    for (size_t l = ix; l < pSim->n(); l += sx)
    {
        double tmp0, tmp1;

        for (size_t i = 0; i < pSim->cstll().log2M(); i++)
        {
            tmp0 = 0.0;
            tmp1 = 0.0;
            for (size_t j = 0; j < pSim->cstll().M(); j++)
            {
                if (pSim->labels()[j] & (1 << (pSim->cstll().log2M() - 1 - i)))
                {
                    tmp1 += exp(-(pSim->mY[l] - pSim->cstll().X()[j]) * (pSim->mY[l] - pSim->cstll().X()[j]) / (2 * pSigma2)) * pSim->cstll().pX()[j];
                }
                else
                {
                    tmp0 += exp(-(pSim->mY[l] - pSim->cstll().X()[j]) * (pSim->mY[l] - pSim->cstll().X()[j]) / (2 * pSigma2)) * pSim->cstll().pX()[j];
                }
            }
            double val = log(tmp0 / tmp1);
            // check usually required when PAS is used with large constellations
            // and severely shaped distributions
            if (isinf(val) == +1)
            {
                llr_tmp[i] = MAX_LLR;
            }
            else if (isinf(val) == -1)
            {
                llr_tmp[i] = MIN_LLR;
            }
            else
            {
                llr_tmp[i] = val;
            }
        }

        for (size_t k = 0; k < pSim->bits(); k++)
        {
            pSim->mLdpcDecoder->mLLRIn[pSim->bit_mapper()[k][l]] = llr_tmp[k];
        }
    }
}

__global__ void cudakernel::sim::calc_llrin(ldpc_sim_device *pSim)
{
    const size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t sx = blockDim.x * gridDim.x;

    for (size_t i = ix; i < pSim->ldpc_code()->nc(); i += sx)
    {
        pSim->mLdpcDecoder->mLLRIn[i] *= (1 - 2 * pSim->mC[i]);
    }
}

__global__ void cudakernel::sim::map_c_to_x(ldpc_sim_device *pSim)
{
    const size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t sx = blockDim.x * gridDim.x;

    size_t tmp;

    for (size_t i = ix; i < pSim->n(); i += sx)
    {
        tmp = 0;
        for (size_t j = 0; j < pSim->bits(); j++)
        {
            tmp += pSim->mC[pSim->bit_mapper()[j][i]] << (pSim->bits() - 1 - j);
        }

        pSim->mX[i] = pSim->labels_rev()[tmp];
    }
}

/*
 *	Decoder kernels
 */
__global__ void cudakernel::decoder::clean_decoder(ldpc_decoder_device *pDecMgd)
{
    const size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t sx = blockDim.x * gridDim.x;

    for (size_t i = ix; i < pDecMgd->mLSum.size(); i += sx)
    {
        pDecMgd->mLSum[i] = 0.0;
        for (size_t l = 0; l < pDecMgd->mLdpcCode->nl(); ++l)
        {
            pDecMgd->mLc2v[l * pDecMgd->mLdpcCode->nnz() + i] = 0.0;
            pDecMgd->mLv2c[l * pDecMgd->mLdpcCode->nnz() + i] = 0.0;
            pDecMgd->mLc2vPre[l * pDecMgd->mLdpcCode->nnz() + i] = 0.0;
        }
    }
}

__global__ void cudakernel::decoder::decode_layered(ldpc_decoder_device *pDecMgd)
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
            pI = pDecMgd->mLdpcCode->nnz() * l;

            //launching kernels
            cudakernel::decoder::decode_lyr_vnupdate<<<gridSizeNC, NUM_THREADS>>>(pDecMgd, pI);
            cudakernel::decoder::decode_lyr_cnupdate<<<get_num_size(pDecMgd->mLdpcCode->lw()[l], NUM_THREADS / 2), NUM_THREADS / 2>>>(pDecMgd, pI, l);
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

__global__ void cudakernel::decoder::decode_lyr_vnupdate(ldpc_decoder_device *pDecMgd, size_t pI)
{
    size_t *vn;
    size_t vw;

    const size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t sx = blockDim.x * gridDim.x;

    //VN processing
    for (size_t i = ix; i < pDecMgd->mLdpcCode->nc(); i += sx)
    {
        double tmp = pDecMgd->mLLRIn[i];
        vw = pDecMgd->mLdpcCode->vn()[i].size();
        vn = pDecMgd->mLdpcCode->vn()[i].data();
        while (vw--)
            tmp += pDecMgd->mLSum[*vn++];

        vn = pDecMgd->mLdpcCode->vn()[i].data();
        vw = pDecMgd->mLdpcCode->vn()[i].size();
        while (vw--)
        {
            pDecMgd->mLv2c[pI + *vn] = tmp - pDecMgd->mLc2v[pI + *vn];
            ++vn;
        }
    }
}

__global__ void cudakernel::decoder::decode_lyr_cnupdate(ldpc_decoder_device *pDecMgd, size_t pI, uint64_t pL)
{
    size_t *cn;
    size_t cw;

    const size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t sx = blockDim.x * gridDim.x;

    double f_tmp[DEC_MAX_DC];
    double b_tmp[DEC_MAX_DC];

    //CN processing
    for (size_t i = ix; i < pDecMgd->mLdpcCode->layers()[pL].size(); i += sx)
    {
        cw = pDecMgd->mLdpcCode->cn()[pDecMgd->mLdpcCode->layers()[pL][i]].size();
        cn = pDecMgd->mLdpcCode->cn()[pDecMgd->mLdpcCode->layers()[pL][i]].data();
        f_tmp[0] = pDecMgd->mLv2c[pI + *cn];
        b_tmp[cw - 1] = pDecMgd->mLv2c[pI + *(cn + cw - 1)];
        for (size_t j = 1; j < cw; j++)
        {
            f_tmp[j] = jacobian(f_tmp[j - 1], pDecMgd->mLv2c[pI + *(cn + j)]);
            b_tmp[cw - 1 - j] = jacobian(b_tmp[cw - j], pDecMgd->mLv2c[pI + *(cn + cw - j - 1)]);
        }

        pDecMgd->mLc2v[pI + *cn] = b_tmp[1];
        pDecMgd->mLc2v[pI + *(cn + cw - 1)] = f_tmp[cw - 2];

        for (size_t j = 1; j < cw - 1; j++)
        {
            pDecMgd->mLc2v[pI + *(cn + j)] = jacobian(f_tmp[j - 1], b_tmp[j + 1]);
        }
    }
}

__global__ void cudakernel::decoder::decode_lyr_sumllr(ldpc_decoder_device *pDecMgd, size_t pI)
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

__global__ void cudakernel::decoder::decode_lyr_appcalc(ldpc_decoder_device *pDecMgd)
{
    size_t *vn;
    size_t vw;

    const size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t sx = blockDim.x * gridDim.x;

    //app calc
    for (size_t i = ix; i < pDecMgd->mLdpcCode->nc(); i += sx)
    {
        pDecMgd->mLLROut[i] = pDecMgd->mLLRIn[i];
        vn = pDecMgd->mLdpcCode->vn()[i].data();
        vw = pDecMgd->mLdpcCode->vn()[i].size();
        while (vw--)
            pDecMgd->mLLROut[i] += pDecMgd->mLSum[*vn++];
        pDecMgd->mCO[i] = (pDecMgd->mLLROut[i] <= 0);
    }
}

__global__ void cudakernel::decoder::calc_synd(ldpc_decoder_device *pDecMgd)
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