#include "ldpcsim.h"

namespace ldpc
{
//Init constructor
constellation::constellation(labels_t pM)
    : mM(pM), mLog2M(log2(pM)), mPX(pM, 1.0 / pM), mX(pM, 0.0)
{
    double m = 0;
    for (std::size_t j = 0; j < mM; ++j)
    {
        mX[j] = (double)-mM + 1 + 2 * j;
        m += mX[j] * mX[j] * mPX[j];
    }

    for (std::size_t j = 0; j < mM; ++j)
    {
        mX[j] = mX[j] / sqrt(m);
    }
}

/*
 * ldpc_sim
 */
//init constructor
ldpc_sim::ldpc_sim(ldpc_code *pCode, const char *pSimFileName, const char *pMapFileName, std::uint16_t numThreads, std::size_t seed)
    : mLdpcCode(pCode), mLdpcDecoder(numThreads, ldpc_decoder(pCode, this, 0, true)), mThreads(numThreads),
    mRNG(numThreads), mRandNormal(numThreads), mRNGSeed(seed)
{
    try
    {
        FILE *fpSim;
        //FILE *fpMap;

        std::ifstream fsSim(pSimFileName);
        //std::ifstream fsMap(pMapFileName);
        std::string fsLine;
        std::string fsSubStr;

        char tmp[1000];
        double m;
        char *tokp;

        std::getline(fsSim, fsLine);
        mLogfile = fsLine.substr(fsLine.find(":") + 2); //name

        labels_t M;
        std::getline(fsSim, fsLine);
        fsSubStr = fsLine.substr(fsLine.find(":") + 2); //M
        M = static_cast<labels_t>(std::stoul(fsSubStr));

        //setup constellation
        mConstellation = constellation(M);

        std::getline(fsSim, fsLine);
        fsSubStr = fsLine.substr(fsLine.find(":") + 2); //bits
        mBits = static_cast<labels_t>(std::stoul(fsSubStr));

        if (mBits != 1 || M != 2)
        {
            throw std::runtime_error("error parsing simfile: the simulation currently only supports BPSK");
        }

        mLabels = std::vector<labels_t>();
        std::getline(fsSim, fsLine);
        fsSubStr = fsLine.substr(fsLine.find(":") + 2); //labels
        strcpy(tmp, fsSubStr.c_str());
        tokp = strtok(tmp, ", ");
        std::size_t i = 0;
        while (tokp != nullptr)
        {
            mLabels.push_back(static_cast<labels_t>(atoi(tokp)));
            tokp = strtok(nullptr, ", ");
            i++;
        }
        if (i != M)
        {
            throw std::runtime_error("error parsing simfile: number of constellation points does not match label size");
        }

        mSnrs = std::vector<double>();
        std::getline(fsSim, fsLine);
        fsSubStr = fsLine.substr(fsLine.find(":") + 2); //snrs
        strcpy(tmp, fsSubStr.c_str());
        tokp = strtok(tmp, ",");
        i = 0;
        while (tokp != nullptr)
        {
            mSnrs.push_back(static_cast<double>(atof(tokp)));
            tokp = strtok(nullptr, ",");
            i++;
        }

        std::getline(fsSim, fsLine);
        fsSubStr = fsLine.substr(fsLine.find(":") + 2); //max frames
        mMaxFrames = std::stoul(fsSubStr);

        std::getline(fsSim, fsLine);
        fsSubStr = fsLine.substr(fsLine.find(":") + 2); //min fec
        mMinFec = std::stoul(fsSubStr);

        std::getline(fsSim, fsLine);
        fsSubStr = fsLine.substr(fsLine.find(":") + 2); //bp iter
        mBPIter = std::stoul(fsSubStr);

        if (mLdpcCode->nct() % mBits != 0)
        {
            throw std::runtime_error("Chosen setting m with n_c does not work. Please correct.");
        }

        mN = mLdpcCode->nct() / mBits;
        mSE = (((double)mLdpcCode->kct()) / mLdpcCode->nct()) * mBits;

        //setup bitmapper
        /*
        mBitMapper = std::vector<std::vector<std::size_t>>(mBits, std::vector<std::size_t>(mN, 0));
        std::getline(fsMap, fsLine);
        std::size_t pos = 0;
        for (auto &ibm : mBitMapper)
        {
            for (auto &jbm : ibm)
            {
                pos = fsLine.find(",");
                fsSubStr = fsLine.substr(0, pos);
                jbm = std::stoul(fsSubStr);
                fsLine.erase(0, pos + 2); //+2 for comma & space
            }
        }
        */
        
        // position of transmitted bits
        mBitPos = std::vector<std::size_t>(mLdpcCode->nct(), 0);
        bool found_p = false;
        bool found_s = false;

        std::size_t idx = 0;
        for (std::size_t i = 0; i < mLdpcCode->nc(); i++)
        {
            for (auto s : mLdpcCode->shorten())
            {
                if (s == i)
                {
                    found_s = true;
                    break;
                }
            }
            for (auto p : mLdpcCode->puncture())
            {
                if (p == i)
                {
                    found_p = true;
                    break;
                }
            }
            if ((!found_s) && (!found_p))
            {
                mBitPos[idx] = i;
                idx++;
            }
            else
            {
                found_p = false;
                found_s = false;
            }
        }

        //channel i/o
        //for many threads we need independent vectors
        mX = mat_double_t(mThreads, vec_double_t(mN, 1.0)); // all one for bpsk and all zero codeword 
        mY = mat_double_t(mThreads, vec_double_t(mN));
        mC = mat_bits_t(mThreads, vec_bits_t(mLdpcCode->nc(), 0)); //all zero

        // RNG setup
        //std::random_device rd;
        //results may vary with same seed, since some threads are executed more than others
        for (std::size_t i = 0; i < mThreads; ++i)
        {
            //decoder
            mLdpcDecoder[i] = ldpc_decoder(mLdpcCode, this, mBPIter, true);

            mRNG[i] = std::mt19937(mRNGSeed + i); // different seeds for threads
            mRandNormal[i] = std::normal_distribution<double>(0.0, 1.0);
        }
    }
    catch (std::exception &e)
    {
        std::cout << "Error: ldpc_sim::ldpc_sim() " << e.what() << "\n";
        exit(EXIT_FAILURE);
    }
}


void ldpc_sim::simulate_awgn(double pSigma2, std::uint16_t threadid)
{
    //double a = 0;
    //double Pn = 0;
    //double Px = 0;

    for (std::size_t i = 0; i < mN; i++)
    {
        //Pn += a * a;
        //Px += mConstellation.X()[mX[i]] * mConstellation.X()[mX[i]];
        mY[threadid][i] = mRandNormal[threadid](mRNG[threadid]) * sqrt(pSigma2) + mX[threadid][i];
    }
}

void ldpc_sim::encode()
{
    /*
    //encode all zero
    for (std::size_t i = 0; i < mLdpcDecoder->nc(); ++i)
    {
        mC[i] = 0;
    }
    */
}

void ldpc_sim::map()
{
    //for higher order modulation we require a bitmapper
    /*
    //for bpsk we have the mapping x_i = 1 - 2*c_i
    for (std::size_t i = 0; i < mN; ++i)
    {
        mX[i] = 1 - 2*mC[mBitPos[i]];
    }
    */
}



void ldpc_sim::print()
{
    std::cout << "logfile: " << mLogfile << "\n";
    printf("n: %lu\n", mN);
    printf("M: %hu\n", mConstellation.M());
    printf("bits: %hu\n", mBits);
    printf("SNRs: ");
    for (const auto &i : mSnrs)
    {
        printf("%lf ", i);
    }
    printf("\n");
    printf("labels:\n");
    for (std::size_t i = 0; i < mConstellation.M(); i++)
    {
        printf("\t%lu: ", i);
        dec2bin(mLabels[i], mBits);
        printf("\n");
    }
    printf("max frames: %lu\n", mMaxFrames);
    printf("min fec: %lu\n", mMinFec);
    printf("bp iter: %lu\n", mLdpcDecoder[0].max_iter());
    printf("early term: %u\n", mLdpcDecoder[0].early_termination());
    printf("SE: %.4lf\n", mSE);
    printf("RNG: mt19937\n");
    printf(" Thread ID | Seed\n");
    for (std::size_t i = 0; i < mThreads; i++)
    {
        printf(" %3d       | %d\n",i, mRNGSeed + i);
    }
}

void ldpc_sim::print_file_header(const char *binaryFile, const char *codeFile, const char *simFile, const char *mapFile)
{
    /*
    FILE *fp = fopen(mLogfile, "a+");
    fprintf(fp, "%% binary: %s (Version: %s, Built: %s)\n", binaryFile, VERSION, BUILD_DATE);
    fprintf(fp, "%% sim file: %s\n", simFile);
    fprintf(fp, "%% code file: %s\n", codeFile);
    fprintf(fp, "%% mapping file: %s\n", mapFile);
    fprintf(fp, "%% result file: %s\n", mLogfile);
    fprintf(fp, "%% iter: %lu\n", mBPIter);
    fprintf(fp, "%% max frames: %lu\n", mMaxFrames);
    fprintf(fp, "%% min fec: %lu\n", mMinFec);
    fprintf(fp, "%% BP early terminate: %hu\n", 1);
    fprintf(fp, "%% num threads: %d\n", 1);
    */
}

void ldpc_sim::log_error(std::size_t pFrameNum, double pSNR)
{
    /*
    char errors_file[MAX_FILENAME_LEN];
    snprintf(errors_file, MAX_FILENAME_LEN, "errors_%s", mLogfile.c_str());

    FILE *fp = fopen(errors_file, "a+");
    if (!fp)


    {


        printf("can not

 open error log file.\n");
        exit(EXIT_FAILU

RE);
    }

    // calculation of syndrome and failed syndrome checks
    std::size_t synd_weight = 0;
    for (auto si : mLdpcDecoder->syndrome())
    {
        synd_weight += si;
    }
    std::vector<std::size_t> failed_checks_idx(synd_weight);
    std::size_t j = 0;
    for (std::size_t i = 0; i < mLdpcCode->mc(); i++)
    {
        if (mLdpcDecoder->syndrome()[i] == 1)
        {
            failed_checks_idx[j++] = i;
        }
    }

    // calculation of failed codeword bits
    std::size_t cw_dis = 0;
    for (std::size_t i = 0; i < mLdpcCode->nc(); i++)
    {
#ifdef ENCODE
        cw_dis += ((mLdpcDecoder->llr_out()[i] <= 0) != mC[pThreads][i]);
#else
        cw_dis += ((mLdpcDecoder->llr_out()[i] <= 0) != 0);
#endif
    }

    std::vector<std::size_t> x(mN);
    std::vector<std::size_t> xhat(mN);
    std::vector<std::size_t> chat(mLdpcCode->nc());

    for (std::size_t i = 0; i < mLdpcCode->nc(); ++i)
    {
        chat[i] = (mLdpcDecoder->llr_out()[i] <= 0);
    }

    std::size_t tmp;
    //map c to x map_c_to_x(c, x);
    for (std::size_t i = 0; i < mN; i++)
    {
        tmp = 0;
        for (std::size_t j = 0; j < mBits; j++)
        {
            tmp += mC[mBitMapper[j][i]] << (mBits - 1 - j);
        }

        x[i] = mLabelsRev[tmp];
    }

    //map_c_to_x(chat, xhat);
    for (std::size_t i = 0; i < mN; i++)
    {
        tmp = 0;
        for (std::size_t j = 0; j < mBits; j++)
        {
            tmp += chat[mBitMapper[j][i]] << (mBits - 1 - j);
        }

        xhat[i] = mLabelsRev[tmp];
    }

    double cw_dis_euc = 0;
    for (std::size_t i = 0; i < mN; i++)
    {
#ifdef ENCODE
        cw_dis_euc += (mConstellation.X()[x[i]] - mConstellation.X()[xhat[i]]) * (mConstellation.X()[x[i]] - mConstellation.X()[xhat[i]]);
#else
        cw_dis_euc += (mConstellation.X()[0] - mConstellation.X()[xhat[i]]) * (mConstellation.X()[0] - mConstellation.X()[xhat[i]]);
#endif
    }
    std::vector<std::size_t> failed_bits_idx(cw_dis);
    j = 0;
    for (std::size_t i = 0; i < mLdpcCode->nc(); i++)
    {
#ifdef ENCODE
        if (chat[i] != mC[pThreads][i])
        {
            failed_bits_idx[j++] = i;
        }
#else
        if (chat[i] != 0)
        {
            failed_bits_idx[j++] = i;
        }
#endif
    }

    // print results in file
    fprintf(fp, "SNR: %.2f -- frame: %lu -- is codeword: %d -- dE(c,chat): %.3f -- dH(c,chat): %lu | ", pSNR, pFrameNum, synd_weight == 0, cw_dis_euc, cw_dis);
    for (auto failed_bits_idx_i : failed_bits_idx)
    {
        fprintf(fp, "%lu ", failed_bits_idx_i);
    }
    fprintf(fp, " -- ");
    fprintf(fp, "synd weight: %lu | ", synd_weight);
    for (auto failed_checks_idx_i : failed_checks_idx)
    {
        fprintf(fp, "%lu ", failed_checks_idx_i);
    }
    fprintf(fp, "\n");
    fclose(fp);
    */
}
} // namespace ldpc
