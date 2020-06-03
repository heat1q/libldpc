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

        mCheckNodeN = mat_size_t(mM);
        mVarNodeN = mat_size_t(mN);
        for (auto &row : mCheckNodeN)
        {
            row = vec_size_t();
        }
        for (auto &col : mVarNodeN)
        {
            col = vec_size_t();
        }

        for (std::size_t i = 0; i < mNNZ; i++)
        {
            fscanf(fpCode, "%lu %lu\n", &(mEdgeCN[i]), &(mEdgeVN[i]));

            // if mEdgeCN[i] is in Shorten skipt both, i.e. the weight of this CN is 0
            // if mEdgeVN[i] is in Puncture or Shorten skip both, i.e. the weight of this VN is 0
            //if (std::find(mShorten.begin(), mShorten.end(), mEdgeCN[i]) == mShorten.end()
            //    && std::find(mPuncture.begin(), mPuncture.end(), mEdgeVN[i]) == mPuncture.end()
            //    && std::find(mShorten.begin(), mShorten.end(), mEdgeVN[i]) == mShorten.end())
            cw[mEdgeCN[i]]++;
            vw[mEdgeVN[i]]++;

            mCheckNodeN[mEdgeCN[i]].push_back(mEdgeVN[i]);
            mVarNodeN[mEdgeVN[i]].push_back(mEdgeCN[i]);
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
    for (auto x : mPuncture)
    {
        std::cout << x << " ";
    }
    std::cout << "\nnum shorten: " << mShorten.size() << "\n";
    std::cout << "shorten: ";
    for (auto x : mShorten)
    {
        std::cout << x << " ";
    }
    std::cout << "\n";
}

std::size_t ldpc_code::calc_rank()
{
    std::size_t rank = mN;
    mat_size_t checkNodeN = mCheckNodeN;
    mat_size_t varNodeN = mVarNodeN;

    for (std::size_t row = 0; row < rank; ++row)
    {
        //std::cout << "Row value: " << row << "\n";

        // check what value h[row][row] has
        auto it = std::find(varNodeN[row].begin(), varNodeN[row].end(), row);
        if (it != varNodeN[row].end()) // values is non-zero
        {
            // now add current row to all rows where a non-zero entry is in the current col, to remove 1
            vec_size_t tmp = varNodeN[row];
            for (std::size_t j = 0; j < tmp.size(); ++j)
            {
                //std::cout << "Check: " << tmp[j] << "\n";
                if (tmp[j] > row)
                {
                    //std::cout << "Add rows " << row << " to " << tmp[j] << "\n";
                    ldpc_code::add_rows(checkNodeN, varNodeN, tmp[j], checkNodeN[row]);
                }
            }
        }
        else // value is zero
        {
            // if there is a row below it with non-zero entry in same col, swap current rows
            bool isZero = true;
            // find first row with non-zero entry
            for (std::size_t j = 0; j < varNodeN[row].size(); ++j)
            {
                if (varNodeN[row][j] > row)
                {
                    //std::cout << "Swap rows " << varNodeN[row][j] << " with " << row << "\n";
                    ldpc_code::swap_rows(checkNodeN, varNodeN, varNodeN[row][j], row);
                    isZero = false;
                    break;
                }
            }

            // if all elements in current col below h[row][row] are zero, swap col it with rank-1 col
            if (isZero)
            {
                --rank;
                // copy last col
                ldpc_code::zero_col(checkNodeN, varNodeN, row);
                ldpc_code::add_cols(checkNodeN, varNodeN, row, varNodeN[rank]);
            }

            --row;
        }
        /*
        std::cout << "CN Perspective:\n";
        for (const auto &vn : checkNodeN)
        {
            pgd::vec_size_t row(this->nc());
            for (auto vn_i : vn)
            {
                row[vn_i] = 1;
            }

            for (std::size_t n = 0; n < this->nc(); ++n)
            {
                std::cout << row[n] << " ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
        */
    }

    return rank;
}

void ldpc_code::swap_rows(mat_size_t &checkNodeN, mat_size_t &varNodeN, std::size_t first, std::size_t second)
{
    vec_size_t first_tmp = checkNodeN[first];
    vec_size_t second_tmp = checkNodeN[second];

    ldpc_code::zero_row(checkNodeN, varNodeN, first);
    ldpc_code::zero_row(checkNodeN, varNodeN, second);

    ldpc_code::add_rows(checkNodeN, varNodeN, first, second_tmp);
    ldpc_code::add_rows(checkNodeN, varNodeN, second, first_tmp);
}

void ldpc_code::swap_cols(mat_size_t &checkNodeN, mat_size_t &varNodeN, std::size_t first, std::size_t second)
{
    vec_size_t first_tmp = varNodeN[first];
    vec_size_t second_tmp = varNodeN[second];

    ldpc_code::zero_col(checkNodeN, varNodeN, first);
    ldpc_code::zero_col(checkNodeN, varNodeN, second);

    ldpc_code::add_cols(checkNodeN, varNodeN, first, second_tmp);
    ldpc_code::add_cols(checkNodeN, varNodeN, second, first_tmp);
}

void ldpc_code::add_rows(mat_size_t &checkNodeN, mat_size_t &varNodeN, std::size_t dest, const vec_size_t &src)
{
    vec_size_t new_row = checkNodeN[dest];
    for (auto vn : src) // append new vn and check if already in
    {
        auto it = std::find(new_row.begin(), new_row.end(), vn);
        if (it == new_row.end())
        {
            new_row.push_back(vn);
        }
        else
        {
            new_row.erase(it);
        }
    }

    ldpc_code::zero_row(checkNodeN, varNodeN, dest); // set row zero

    checkNodeN[dest] = new_row;

    // append to vn
    for (auto vn : new_row)
    {
        varNodeN[vn].push_back(dest);
    }
}

void ldpc_code::add_cols(mat_size_t &checkNodeN, mat_size_t &varNodeN, std::size_t dest, const vec_size_t &src)
{
    vec_size_t new_col = varNodeN[dest];
    for (auto cn : src) // append new cn and check if already in
    {
        auto it = std::find(new_col.begin(), new_col.end(), cn);
        if (it == new_col.end())
        {
            new_col.push_back(cn);
        }
        else
        {
            new_col.erase(it);
        }
    }

    ldpc_code::zero_col(checkNodeN, varNodeN, dest); // set row zero

    varNodeN[dest] = new_col;

    // append to cn
    for (auto cn : new_col)
    {
        checkNodeN[cn].push_back(dest);
    }
}

void ldpc_code::zero_row(mat_size_t &checkNodeN, mat_size_t &varNodeN, std::size_t m)
{
    for (auto vn : checkNodeN[m]) // from selected row, for each vn index, remove m from vn
    {
        varNodeN[vn].erase(std::remove(varNodeN[vn].begin(), varNodeN[vn].end(), m), varNodeN[vn].end());
    }
    checkNodeN[m] = vec_size_t();
}

void ldpc_code::zero_col(mat_size_t &checkNodeN, mat_size_t &varNodeN, std::size_t n)
{
    for (auto cn : varNodeN[n]) // from selected col, for each cn index, remove n from cn
    {
        checkNodeN[cn].erase(std::remove(checkNodeN[cn].begin(), checkNodeN[cn].end(), n), checkNodeN[cn].end());
    }
    varNodeN[n] = vec_size_t();
}

} // namespace pgd
