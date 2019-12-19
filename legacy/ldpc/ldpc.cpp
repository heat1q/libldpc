#include "ldpc.h"
#include <algorithm>
#include <iterator>

namespace pgd
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
        vec_size_t cw(mM, 0);
        vec_size_t vw(mN, 0);

        mEdgeCN = vec_size_t(mNNZ);
        mEdgeVN = vec_size_t(mNNZ);

        for (std::size_t i = 0; i < mNNZ; i++)
        {
            fscanf(fpCode, "%lu %lu\n", &(mEdgeCN[i]), &(mEdgeVN[i]));

            // if mEdgeCN[i] is in Shorten skipt both, i.e. the weight of this CN is 0
            // if mEdgeVN[i] is in Puncture or Shorten skip both, i.e. the weight of this VN is 0
            //if (std::find(mShorten.begin(), mShorten.end(), mEdgeCN[i]) == mShorten.end()
            //    && std::find(mPuncture.begin(), mPuncture.end(), mEdgeVN[i]) == mPuncture.end() 
            //    && std::find(mShorten.begin(), mShorten.end(), mEdgeVN[i]) == mShorten.end())
            {
                cw[mEdgeCN[i]]++;
                vw[mEdgeVN[i]]++;
            }
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
            mCN[mEdgeCN[i]][cwTmp[mEdgeCN[i]]++] = i;
            mVN[mEdgeVN[i]][vwTmp[mEdgeVN[i]]++] = i;
        }

        // maximum check node degree
        mMaxDC = *(std::max_element(cw.begin(), cw.end()));

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
    std::cout << "puncture: ";
    for (auto x: mPuncture)
    {
        std::cout << x << " ";
    }
    std::cout << "\nnum shorten: " << mShorten.size() << "\n";
    std::cout << "shorten: ";
    for (auto x: mShorten)
    {
        std::cout << x << " ";
    }
    std::cout << "\n";
}
} // namespace pgd
