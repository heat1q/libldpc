#include "../src/sim/ldpcsim.h"

//layered cpu decoder
std::size_t test_decode_layered_cpu(ldpc::cuda_ptr<ldpc::ldpc_decoder> &dec, std::size_t iter)
{
    std::size_t *vn;
    std::size_t *cn;

    std::size_t vw;
    std::size_t cw;

    double mExMsgCN[100];

    //initialize
    for (std::size_t i = 0; i < dec->mLdpcCode->nnz(); ++i)
    {
        dec->mLv2c[i] = dec->mLLRIn[dec->mLdpcCode->c()[i]];
        dec->mLc2v[i] = 0.0;
    }

    std::size_t I = 0;
    while (I < iter)
    {
        for (std::size_t l = 0; l < dec->mLdpcCode->layers().size(); ++l)
        {
            // CN processing
            for (std::size_t i = 0; i < dec->mLdpcCode->layers()[l].size(); ++i)
            {
                cw = dec->mLdpcCode->cn()[dec->mLdpcCode->layers()[l][i]].size();
                cn = dec->mLdpcCode->cn()[dec->mLdpcCode->layers()[l][i]].get();

                double tmp = 1;
                for (std::size_t j = 0; j < cw; ++j)
                {
                    mExMsgCN[j] = 1 - 2 / (exp(dec->mLv2c[cn[j]]) + 1); //tanh(mLv2c[cn[j]]);
                    tmp *= mExMsgCN[j];
                }

                for (std::size_t j = 0; j < cw; ++j)
                {
                    dec->mLc2v[cn[j]] = log((mExMsgCN[j] + tmp) / (mExMsgCN[j] - tmp)); //2*atanh(tmp/mExMsgCN[j]);
                }
            }

            // VN processing and app calc
            for (std::size_t i = 0; i < dec->mLdpcCode->nc(); ++i)
            {
                double tmp = dec->mLLRIn[i];
                vw = dec->mLdpcCode->vn()[i].size(); // degree of VN
                vn = dec->mLdpcCode->vn()[i].get();  //neighbours of VN
                while (vw--)
                    tmp += dec->mLc2v[*vn++];

                dec->mCO[i] = (dec->mLLROut[i] <= 0); // approx decision on ith bits
                dec->mLLROut[i] = tmp;

                vw = dec->mLdpcCode->vn()[i].size(); // degree of VN
                vn = dec->mLdpcCode->vn()[i].get();  //neighbours of VN
                while (vw--)
                {
                    dec->mLv2c[*vn] = tmp - dec->mLc2v[*vn];
                    ++vn;
                }
            }

            if (dec->early_termination())
            {
                printf("Should not go here!");
            }
        }

        ++I;
    }

    return I;
}

__global__ void test_kern_decode(ldpc::cuda_ptr<ldpc::ldpc_decoder> *dec_vec, ldpc::labels_t iter)
{
    ldpc::labels_t ix = threadIdx.x;
    ldpc::ldpc_decoder *dec = dec_vec[ix].get();

    ldpc::cudakernel::decoder::init_decoder<<<ldpc::get_num_size(dec->mLdpcCode->nc(), NUMK_THREADS), NUMK_THREADS>>>(dec);
    ldpc::labels_t I = 0;
    while (I < iter)
    {
        for (std::size_t l = 0; l < dec->mLdpcCode->nl(); ++l)
        {
            //launching kernels
            ldpc::cudakernel::decoder::decode_lyr_cnupdate<<<ldpc::get_num_size(dec->mLdpcCode->lw()[l], NUMK_THREADS), NUMK_THREADS>>>(dec, l);
            ldpc::cudakernel::decoder::decode_lyr_appcalc<<<ldpc::get_num_size(dec->mLdpcCode->nc(), NUMK_THREADS), NUMK_THREADS>>>(dec);

            if (dec->early_termination())
            {
                printf("Should not go here!");
            }
        }

        ++I;
    }
}

__global__ void test_kern_iter_tp(ldpc::cuda_ptr<ldpc::ldpc_decoder> *dec_vec, std::size_t dec_len, int clk, double *time, double *tp)
{

    for (ldpc::labels_t iter = 1; iter <= 200; ++iter)
    {
        time[iter - 1] = 1e12;
        for (std::size_t k = 0; k < 15; ++k)
        {
            double starttime = clock64();
            //decode
            test_kern_decode<<<1, dec_len>>>(dec_vec, iter);
            cudaDeviceSynchronize();

            double endtime = clock64();
            time[iter - 1] = fmin(time[iter - 1], fabs((endtime - starttime) / (clk * dec_len)));
            tp[iter - 1] = 1e-3 * dec_vec[0]->mLdpcCode->nc() / time[iter - 1];
        }
        printf("GPU: Code Length: %u -- Iterations: %lu -- Time: %.3f -- Throughput: %.2f Mbits/s\n", dec_vec[0]->mLdpcCode->nc(), iter, time[iter - 1], tp[iter - 1]);
    }
}

__global__ void test_kern_bl_tp(ldpc::cuda_ptr<ldpc::ldpc_decoder> *dec_vec, std::size_t dec_len, int clk, double *time, double *tp)
{
    *time = 1e12;
    for (std::size_t k = 0; k < 15; ++k)
    {
        double starttime = clock64();
        //decode
        test_kern_decode<<<1, dec_len>>>(dec_vec, dec_vec[0]->max_iter());
        cudaDeviceSynchronize();

        double endtime = clock64();
        *time = fmin(*time, fabs((endtime - starttime) / (clk * dec_len)));
        *tp = 1e-3 * dec_vec[0]->mLdpcCode->nc() / (*time);
    }
    printf("GPU: Code Length: %u -- Iterations: %lu -- Time: %.3f -- Throughput: %.2f Mbits/s\n", dec_vec[0]->mLdpcCode->nc(), dec_vec[0]->max_iter(), *time, *tp);
}

void test_tp_iter(ldpc::cuda_vector<ldpc::cuda_ptr<ldpc::ldpc_decoder>> &dec_vec)
{
    // create LLRs
    for (auto &dec : dec_vec)
    {
        for (auto &x : dec->mLLRIn)
        {
            x = ldpc::ldpc_sim_device::randn() * 2;
        }
    }
    ldpc::cuda_vector<double> time(200);
    ldpc::cuda_vector<double> tp(200);

    int peak_clk = 1;
    cudaDeviceGetAttribute(&peak_clk, cudaDevAttrClockRate, 0);

    test_kern_iter_tp<<<1, 1>>>(dec_vec.get(), dec_vec.size(), peak_clk, time.get(), tp.get());
    cudaDeviceSynchronize();

    std::ofstream fp;
    fp.open("res_tp_iter_gpu.txt");
    fp << "iter time_ms tp_mbits\n";
    for (std::size_t i = 0; i < 200; ++i)
    {
        fp << i + 1 << " " << time[i] << " " << tp[i] << "\n";
    }
    fp.close();

    // test cpu
    for (std::size_t iter = 1; iter <= 200; ++iter)
    {
        time[iter - 1] = 1e12;
        for (std::size_t k = 0; k < 20; ++k)
        {
            auto start = std::chrono::high_resolution_clock::now();
            test_decode_layered_cpu(dec_vec[0], iter);
            auto elapsed = std::chrono::high_resolution_clock::now() - start;
            time[iter - 1] = fmin(time[iter - 1], fabs(static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed).count()) * 1e-6));
            tp[iter - 1] = 1e-3 * dec_vec[0]->mLdpcCode->nc() / time[iter - 1];
        }
        printf("CPU: Code Length: %u -- Iterations: %lu -- Time: %.3f -- Throughput: %.2f Mbits/s\n", dec_vec[0]->mLdpcCode->nc(), iter, time[iter - 1], tp[iter - 1]);
    }

    fp.open("res_tp_iter_cpu.txt");
    fp << "iter time_ms tp_mbits\n";
    for (std::size_t i = 0; i < 200; ++i)
    {
        fp << i + 1 << " " << time[i] << " " << tp[i] << "\n";
    }
    fp.close();
}
/*
tp vs iter: (fixed n)
    argv[1]: codefile
    argv[2]: num frames in parallel
tp vs n: (fixed iter)
    argv[1]: iter
    argv[2]: num samples
    argv[3]: num frames in parallel
 */
int main(int argc, char *argv[])
{
    if (argc == 3)
    {
        ldpc::cuda_ptr<ldpc::ldpc_code_device> code_dev(
            ldpc::ldpc_code_device(
                argv[1], ""));

        ldpc::cuda_vector<ldpc::cuda_ptr<ldpc::ldpc_decoder>> dec_vec;
        for (std::size_t i = 0; i < atoi(argv[2]); i++)
        {
            dec_vec.push_back(ldpc::cuda_ptr<ldpc::ldpc_decoder>(ldpc::ldpc_decoder(code_dev, 10, false)));
        }

        test_tp_iter(dec_vec);
    }
    else if (argc == 4)
    {
        int samples = atoi(argv[2]);

        ldpc::cuda_vector<double> blocklength(samples);
        ldpc::cuda_vector<double> time_gpu(samples);
        ldpc::cuda_vector<double> time_cpu(samples);
        ldpc::cuda_vector<double> tp_gpu(samples);
        ldpc::cuda_vector<double> tp_cpu(samples);
        int peak_clk = 1;
        cudaDeviceGetAttribute(&peak_clk, cudaDevAttrClockRate, 0);
        char str[128];
        std::ofstream fp;

        for (std::size_t i = 0; i < samples; ++i)
        {
            sprintf(str, "dat/code_dv3_dc6_i=%u.txt", i);
            ldpc::cuda_ptr<ldpc::ldpc_code_device> code_dev(ldpc::ldpc_code_device(str, ""));
            ldpc::cuda_vector<ldpc::cuda_ptr<ldpc::ldpc_decoder>> dec_vec;
            for (std::size_t i = 0; i < atoi(argv[3]); i++)
            {
                dec_vec.push_back(ldpc::cuda_ptr<ldpc::ldpc_decoder>(ldpc::ldpc_decoder(code_dev, atoi(argv[1]), false)));
            }

            for (auto &dec : dec_vec)
            {
                for (auto &x : dec->mLLRIn)
                {
                    x = ldpc::ldpc_sim_device::randn() * 2;
                }
            }

            test_kern_bl_tp<<<1, 1>>>(dec_vec.get(), dec_vec.size(), peak_clk, time_gpu.get() + i, tp_gpu.get() + i);
            cudaDeviceSynchronize();

            //test cpu
            time_cpu[i] = 1e12;
            for (std::size_t k = 0; k < 25; ++k)
            {
                auto start = std::chrono::high_resolution_clock::now();
                test_decode_layered_cpu(dec_vec[0], dec_vec[0]->max_iter());
                auto elapsed = std::chrono::high_resolution_clock::now() - start;
                time_cpu[i] = fmin(time_cpu[i], fabs(static_cast<double>(std::chrono::duration_cast<std::chrono::nanoseconds>(elapsed).count()) * 1e-6));
                tp_cpu[i] = 1e-3 * dec_vec[0]->mLdpcCode->nc() / time_cpu[i];
            }
            printf("CPU: Code Length: %u -- Iterations: %lu -- Time: %.3f -- Throughput: %.2f Mbits/s\n", dec_vec[0]->mLdpcCode->nc(), dec_vec[0]->max_iter(), time_cpu[i], tp_cpu[i]);
            blocklength[i] = code_dev->nc();
        }

        fp.open("res_tp_bl.txt");
        fp << "n gpu_time_ms gpu_tp_mbits cpu_time_ms cpu_tp_mbits\n";
        for (std::size_t i = 0; i < samples; ++i)
        {
            fp << blocklength[i] << " " << time_gpu[i] << " " << tp_gpu[i] << " " << time_cpu[i] << " " << tp_cpu[i] << "\n";
        }
        fp.close();
    }
}