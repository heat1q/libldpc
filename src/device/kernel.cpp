#include "../sim/ldpcsim.h"

namespace ldpc
{
/*
*	Simulation kernels
*/
__global__ void cudakernel::sim::setup_rng(ldpc_sim_device *pSim)
{
    // initialize curand
    const labels_t ix = blockIdx.x;
    const std::size_t tx = threadIdx.x;
    const std::size_t sx = blockDim.x;

    for (std::size_t i = tx; i < pSim->n(); i += sx)
    {
        curand_init(pSim->threads()*ix + i, i, 0, &(pSim->mCurandState[ix][i]));
    }

    for (std::size_t j = tx; j < pSim->mLdpcCode->nct(); j += sx)
    {
        curand_init(pSim->threads()*ix + j, j, 0, &(pSim->mCurandStateEncoding[ix][j]));
    }
}

__global__ void cudakernel::sim::frame_proc(ldpc_sim_device *pSim, double pSigma2)
{
    std::size_t pI;
    const labels_t ix = blockIdx.x;
    ldpc_decoder *pDecMgd = pSim->mLdpcDecoderVec[ix].get();
    const std::size_t gridSizeNC = get_num_size(pDecMgd->mLdpcCode->nc(), NUMK_THREADS);
    const std::size_t gridSizeNNZ = get_num_size(pDecMgd->mLdpcCode->nnz(), NUMK_THREADS);

    //encodeall0
    cudakernel::sim::encode_all0<<<get_num_size(pSim->mLdpcCode->nct(), NUMK_THREADS), NUMK_THREADS>>>(pSim, ix);

    //map c to x
    cudakernel::sim::map_c_to_x<<<get_num_size(pSim->n(), NUMK_THREADS), NUMK_THREADS>>>(pSim, ix);

    //sim awgn
    cudakernel::sim::awgn<<<get_num_size(pSim->n(), NUMK_THREADS), NUMK_THREADS>>>(pSim, pSigma2, ix);

    //num & punture
    //TODO

    //calc llrin
    cudakernel::sim::calc_llrs<<<get_num_size(pSim->n(), NUMK_THREADS), NUMK_THREADS>>>(pSim, pSigma2, ix);
    cudakernel::sim::calc_llrin<<<get_num_size(pSim->mLdpcCode->nc(), NUMK_THREADS), NUMK_THREADS>>>(pSim, ix);

    //decode
    //zero everything out
    cudakernel::decoder::init_decoder<<<gridSizeNNZ, NUMK_THREADS>>>(pDecMgd);

    labels_t I = 0;
    while (I < pDecMgd->max_iter())
    {
        for (std::size_t l = 0; l < pDecMgd->mLdpcCode->nl(); ++l)
        {
            //launching kernels
            cudakernel::decoder::decode_lyr_cnupdate<<<get_num_size(pDecMgd->mLdpcCode->lw()[l], NUMK_THREADS), NUMK_THREADS>>>(pDecMgd, l);
            cudakernel::decoder::decode_lyr_appcalc<<<gridSizeNC, NUMK_THREADS>>>(pDecMgd);
            //printf("in: %.3f    out: %.3f\n",pDecMgd->mLLRIn[0],pDecMgd->mLLROut[0]);

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

    pSim->mLdpcDecoderVec[ix]->mIter = I;
}

//measure the constant time for frame processing over pCount samples
#ifdef LOG_TP
__global__ void cudakernel::sim::frame_time(ldpc_sim_device *pSim, double pSigma2)
{
    std::size_t pI;
    const labels_t ix = blockIdx.x;
    ldpc_decoder *pDecMgd = pSim->mLdpcDecoderVec[ix].get();
    const std::size_t gridSizeNC = get_num_size(pDecMgd->mLdpcCode->nc(), NUMK_THREADS);
    const std::size_t gridSizeNNZ = get_num_size(pDecMgd->mLdpcCode->nnz(), NUMK_THREADS);

    cudakernel::sim::encode_all0<<<get_num_size(pSim->mLdpcCode->nct(), NUMK_THREADS), NUMK_THREADS>>>(pSim, ix);
    cudakernel::sim::map_c_to_x<<<get_num_size(pSim->n(), NUMK_THREADS), NUMK_THREADS>>>(pSim, ix);
    cudakernel::sim::awgn<<<get_num_size(pSim->n(), NUMK_THREADS), NUMK_THREADS>>>(pSim, pSigma2, ix);
    cudakernel::sim::calc_llrs<<<get_num_size(pSim->n(), NUMK_THREADS), NUMK_THREADS>>>(pSim, pSigma2, ix);
    cudakernel::sim::calc_llrin<<<get_num_size(pSim->mLdpcCode->nc(), NUMK_THREADS), NUMK_THREADS>>>(pSim, ix);
}
#endif

//encode the input to all zero
__global__ void cudakernel::sim::encode_all0(ldpc_sim_device *pSim, labels_t pBlockID)
{
    const std::size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const std::size_t sx = blockDim.x * gridDim.x;

    for (std::size_t i = ix; i < pSim->mLdpcCode->nct(); i += sx)
    {
        pSim->mC[pBlockID][pSim->bits_pos()[i]] = (curand_uniform(&(pSim->mCurandStateEncoding[pBlockID][i])) > 0.5000);
    }

    for (std::size_t i = ix; i < pSim->mLdpcCode->num_puncture(); i += sx)
    {
        pSim->mC[pBlockID][pSim->mLdpcCode->puncture()[i]] = (curand_uniform(&(pSim->mCurandStateEncoding[pBlockID][i])) > 0.5000);
    }

    for (std::size_t i = 0; i < pSim->mLdpcCode->num_shorten(); i++)
    {
        pSim->mC[pBlockID][pSim->mLdpcCode->shorten()[i]] = 0;
    }

    //map_c_to_x()
}

__global__ void cudakernel::sim::awgn(ldpc_sim_device *pSim, double pSigma2, labels_t pBlockID)
{
    const std::size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const std::size_t sx = blockDim.x * gridDim.x;

    double a = 0;

    for (std::size_t i = ix; i < pSim->n(); i += sx)
    {
        a = curand_normal(&(pSim->mCurandState[pBlockID][i])) * sqrt(pSigma2);
        pSim->mY[pBlockID][i] = pSim->cstll().X()[pSim->mX[pBlockID][i]] + a;
        //Pn += a * a;
        //Px += pSim->cstll().X()[pSim->mX[pBlockID][i]] * pSim->cstll().X()[pSim->mX[pBlockID][i]];
    }

    //return Px/Pn
}

__global__ void cudakernel::sim::calc_llrs(ldpc_sim_device *pSim, double pSigma2, labels_t pBlockID)
{
    const std::size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const std::size_t sx = blockDim.x * gridDim.x;

    double llr_tmp[SIM_NUM_BITS];

    for (std::size_t l = ix; l < pSim->n(); l += sx)
    {
        double tmp0, tmp1;

        for (std::size_t i = 0; i < pSim->cstll().log2M(); i++)
        {
            tmp0 = 0.0;
            tmp1 = 0.0;
            for (std::size_t j = 0; j < pSim->cstll().M(); j++)
            {
                if (pSim->labels()[j] & (1 << (pSim->cstll().log2M() - 1 - i)))
                {
                    tmp1 += exp(-(pSim->mY[pBlockID][l] - pSim->cstll().X()[j]) * (pSim->mY[pBlockID][l] - pSim->cstll().X()[j]) / (2 * pSigma2)) * pSim->cstll().pX()[j];
                }
                else
                {
                    tmp0 += exp(-(pSim->mY[pBlockID][l] - pSim->cstll().X()[j]) * (pSim->mY[pBlockID][l] - pSim->cstll().X()[j]) / (2 * pSigma2)) * pSim->cstll().pX()[j];
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

        for (std::size_t k = 0; k < pSim->bits(); k++)
        {
            pSim->mLdpcDecoderVec[pBlockID]->mLLRIn[pSim->bit_mapper()[k][l]] = llr_tmp[k];
        }
    }
}

__global__ void cudakernel::sim::calc_llrin(ldpc_sim_device *pSim, labels_t pBlockID)
{
    const std::size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const std::size_t sx = blockDim.x * gridDim.x;

    for (std::size_t i = ix; i < pSim->mLdpcCode->nc(); i += sx)
    {
        pSim->mLdpcDecoderVec[pBlockID]->mLLRIn[i] *= (1 - 2 * pSim->mC[pBlockID][i]);
    }
}

__global__ void cudakernel::sim::map_c_to_x(ldpc_sim_device *pSim, labels_t pBlockID)
{
    const std::size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const std::size_t sx = blockDim.x * gridDim.x;

    std::size_t tmp;

    for (std::size_t i = ix; i < pSim->n(); i += sx)
    {
        tmp = 0;
        for (std::size_t j = 0; j < pSim->bits(); j++)
        {
            tmp += pSim->mC[pBlockID][pSim->bit_mapper()[j][i]] << (pSim->bits() - 1 - j);
        }

        pSim->mX[pBlockID][i] = pSim->labels_rev()[tmp];
    }
}

/*
 *	Decoder kernels
 */
__global__ void cudakernel::decoder::init_decoder(ldpc_decoder *pDecMgd)
{
    const std::size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const std::size_t sx = blockDim.x * gridDim.x;

    for (std::size_t i = ix; i < pDecMgd->mLdpcCode->nnz(); i += sx)
    {
        pDecMgd->mLc2v[i] = 0.0;
        pDecMgd->mLv2c[i] = pDecMgd->mLLRIn[pDecMgd->mLdpcCode->c()[i]];
    }
}

__global__ void cudakernel::decoder::decode_lyr_cnupdate(ldpc_decoder *pDecMgd, std::size_t pL)
{
    std::size_t *cn;
    std::size_t cw;

    const std::size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const std::size_t sx = blockDim.x * gridDim.x;

    double f_tmp[DEC_MAX_DC];
    double b_tmp[DEC_MAX_DC];

    //CN processing
    for (std::size_t i = ix; i < pDecMgd->mLdpcCode->layers()[pL].size(); i += sx)
    {
        cw = pDecMgd->mLdpcCode->cn()[pDecMgd->mLdpcCode->layers()[pL][i]].size();
        cn = pDecMgd->mLdpcCode->cn()[pDecMgd->mLdpcCode->layers()[pL][i]].get();

        cw = pDecMgd->mLdpcCode->cn()[pDecMgd->mLdpcCode->layers()[pL][i]].size();
        cn = pDecMgd->mLdpcCode->cn()[pDecMgd->mLdpcCode->layers()[pL][i]].get();
        f_tmp[0] = pDecMgd->mLv2c[*cn];
        b_tmp[cw - 1] = pDecMgd->mLv2c[*(cn + cw - 1)];
        for (std::size_t j = 1; j < cw; j++)
        {
            f_tmp[j] = jacobian(f_tmp[j - 1], pDecMgd->mLv2c[*(cn + j)]);
            b_tmp[cw - 1 - j] = jacobian(b_tmp[cw - j], pDecMgd->mLv2c[*(cn + cw - j - 1)]);
        }

        pDecMgd->mLc2v[*cn] = b_tmp[1];
        pDecMgd->mLc2v[*(cn + cw - 1)] = f_tmp[cw - 2];

        for (std::size_t j = 1; j < cw - 1; j++)
        {
            pDecMgd->mLc2v[*(cn + j)] = jacobian(f_tmp[j - 1], b_tmp[j + 1]);
        }
    }
}

__global__ void cudakernel::decoder::decode_lyr_appcalc(ldpc_decoder *pDecMgd)
{
    std::size_t *vn;
    std::size_t vw;

    const std::size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const std::size_t sx = blockDim.x * gridDim.x;

    //app calc
    for (std::size_t i = ix; i < pDecMgd->mLdpcCode->nc(); i += sx)
    {
        pDecMgd->mLLROut[i] = pDecMgd->mLLRIn[i];
        vn = pDecMgd->mLdpcCode->vn()[i].get();  //neighbours of VN
        vw = pDecMgd->mLdpcCode->vn()[i].size(); // degree of VN
        while (vw--)
            pDecMgd->mLLROut[i] += pDecMgd->mLc2v[*vn++];

        pDecMgd->mCO[i] = (pDecMgd->mLLROut[i] <= 0); // approx decision on ith bits

        vn = pDecMgd->mLdpcCode->vn()[i].get();  //neighbours of VN
        vw = pDecMgd->mLdpcCode->vn()[i].size(); // degree of VN
        while (vw--)
        {
            pDecMgd->mLv2c[*vn] = pDecMgd->mLLROut[i] - pDecMgd->mLc2v[*vn];
            ++vn;
        }
    }
}

__global__ void cudakernel::decoder::calc_synd(ldpc_decoder *pDecMgd)
{
    const std::size_t ix = blockIdx.x * blockDim.x + threadIdx.x;
    const std::size_t sx = blockDim.x * gridDim.x;

    for (std::size_t i = ix; i < pDecMgd->mLdpcCode->mc(); i += sx)
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
} // namespace ldpc
