#include "ldpc.h"

namespace ldpc
{
/**
 * @brief Construct a new ldpc code::ldpc code object
 * 
 * @param pFileName 
 */
ldpc_code::ldpc_code(const char *pFileName)
    : mMaxDC(0)
{
    try
    {
        FILE *fpCode = fopen(pFileName, "r");
        if (!fpCode)
        {
            throw std::runtime_error("can not open codefile for reading.");
        }

        fscanf(fpCode, "nc: %lu\n", &mN);
        fscanf(fpCode, "mc: %lu\n", &mM);
        fscanf(fpCode, "nct: %lu\n", &mNCT);
        fscanf(fpCode, "mct: %lu\n", &mMCT);
        fscanf(fpCode, "nnz: %lu\n", &mNNZ);
        mK = mN - mM;
        mKCT = mNCT - mMCT;

        std::size_t numPuncture = 0;
        std::size_t numShorten = 0;

        fscanf(fpCode, "puncture [%lu]: ", &numPuncture);
        if (numPuncture != 0)
        {
            mPuncture = vec_size_t(numPuncture);
            for (std::size_t i = 0; i < numPuncture; i++)
            {
                fscanf(fpCode, " %lu ", &(mPuncture[i]));
            }
        }

        fscanf(fpCode, "shorten [%lu]: ", &numShorten);
        if (numShorten != 0)
        {
            mShorten = vec_size_t(numShorten);
            for (std::size_t i = 0; i < numShorten; i++)
            {
                fscanf(fpCode, " %lu ", &(mShorten[i]));
            }
        }

        vec_size_t cwTmp(mM);
        vec_size_t vwTmp(mN);
        vec_size_t cw(mM);
        vec_size_t vw(mN);

        mR = vec_size_t(mNNZ);
        mC = vec_size_t(mNNZ);

        for (std::size_t i = 0; i < mNNZ; i++)
        {
            fscanf(fpCode, "%lu %lu\n", &(mR[i]), &(mC[i]));
            cw[mR[i]]++;
            vw[mC[i]]++;
        }

        mCN = mat_size_t(mM, vec_size_t());
        for (std::size_t i = 0; i < mM; i++)
        {
            mCN[i] = vec_size_t(cw[i]);
        }

        mVN = mat_size_t(mN, vec_size_t());
        for (std::size_t i = 0; i < mN; i++)
        {
            mVN[i] = vec_size_t(vw[i]);
        }

        for (std::size_t i = 0; i < mNNZ; i++)
        {
            mCN[mR[i]][cwTmp[mR[i]]++] = i;
            mVN[mC[i]][vwTmp[mC[i]]++] = i;
        }

        for (std::size_t i = 0; i < mM; i++)
        {
            if (cw[i] > mMaxDC)
            {
                mMaxDC = cw[i];
            }
        }

        fclose(fpCode);
    }
    catch (std::exception &e)
    {
        std::cout << "Error: ldpc_code(): " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
}

/**
 * @brief Prints parameters of LDPC code
 * 
 */
void ldpc_code::print()
{
    std::cout << "nc : " << mN << "\n";
    std::cout << "mc : " << mM << "\n";
    std::cout << "kc : " << mK << "\n";
    std::cout << "nnz : " << mNNZ << "\n";
    std::cout << "nct :" << mNCT << "\n";
    std::cout << "mct : " << mMCT << "\n";
    std::cout << "kct : " << mKCT << "\n";
    std::cout << "max dc : " << mMaxDC << "\n";
    std::cout << "num puncture: " << mPuncture.size() << "\n";
    std::cout << "num shorten: " << mShorten.size() << "\n";
}

/**
 * @brief Prints a decimal as binary with m bits
 * 
 * @param val 
 * @param m 
 */
void dec2bin(std::size_t val, uint8_t m)
{
    for (std::size_t i = 0; i < m; i++)
    {
        printf("%lu", (val >> (m - i - 1) & 0x01));
    }
}

/**
 * @brief Sign function
 * 
 * @param a 
 * @return int 
 */
int sign(double a)
{
    return (a <= 0) ? -1 : 1;
}
} // namespace ldpc