#include "ldpc.h"
#include <iterator>

namespace ldpc
{
/**
 * @brief Construct a new ldpc code::ldpc code object
 * 
 * @param pFileName 
 */
ldpc_code::ldpc_code(const std::string &pFileName)
    : mMaxDC(0)
{
    try
    {
        FILE *fpCode = fopen(pFileName.c_str(), "r");
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

        u64 numPuncture = 0;
        u64 numShorten = 0;

        fscanf(fpCode, "puncture [%lu]: ", &numPuncture);
        if (numPuncture != 0)
        {
            mPuncture = vec_u64(numPuncture);
            for (u64 i = 0; i < numPuncture; i++)
            {
                fscanf(fpCode, " %lu ", &(mPuncture[i]));
            }
        }

        fscanf(fpCode, "shorten [%lu]: ", &numShorten);
        if (numShorten != 0)
        {
            mShorten = vec_u64(numShorten);
            for (u64 i = 0; i < numShorten; i++)
            {
                fscanf(fpCode, " %lu ", &(mShorten[i]));
            }
        }

        mEdgeCN = vec_u64(mNNZ);
        mEdgeVN = vec_u64(mNNZ);

        mCN = mat_u64(mM, vec_u64());
        mVN = mat_u64(mN, vec_u64());

        for (u64 i = 0; i < mNNZ; i++)
        {
            // read the non-zero entries, i.e. edges
            fscanf(fpCode, "%lu %lu\n", &(mEdgeCN[i]), &(mEdgeVN[i]));

            // save edge index to coressponding CN & VN
            mCN[mEdgeCN[i]].push_back(i);
            mVN[mEdgeVN[i]].push_back(i);
        }

        // maximum check node degree
        auto tmp = std::max_element(mCN.begin(), mCN.end(), [](const vec_u64 &a, const vec_u64 &b) { return (a.size() < b.size()); });
        mMaxDC = tmp->size();

        // position of transmitted bits
        for (u64 i = 0; i < mN; i++)
        {
            auto tmp = std::find(mShorten.cbegin(), mShorten.cend(), i);
            if (tmp != mShorten.cend()) continue; // skip if current index shortened
            tmp = std::find(mPuncture.cbegin(), mPuncture.cend(), i);
            if (tmp != mPuncture.cend()) continue; // skip if current index punctured

            mBitPos.push_back(i);
        }

        fclose(fpCode);
    }
    catch (std::exception &e)
    {
        std::cout << "Error: ldpc_code(): " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
}


u64 ldpc_code::calc_rank()
{
    u64 rank = mN;
    mat_u64 checkNodeN = mCheckNodeN;
    mat_u64 varNodeN = mVarNodeN;

    for (u64 row = 0; row < rank; ++row)
    {
        //std::cout << "Row value: " << row << "\n";

        // check what value h[row][row] has
        auto it = std::find(varNodeN[row].begin(), varNodeN[row].end(), row);
        if (it != varNodeN[row].end()) // values is non-zero
        {
            // now add current row to all rows where a non-zero entry is in the current col, to remove 1
            vec_u64 tmp = varNodeN[row];
            for (u64 j = 0; j < tmp.size(); ++j)
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
            for (u64 j = 0; j < varNodeN[row].size(); ++j)
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
            ldpc::vec_u64 row(this->nc());
            for (auto vn_i : vn)
            {
                row[vn_i] = 1;
            }

            for (u64 n = 0; n < this->nc(); ++n)
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

void ldpc_code::swap_rows(mat_u64 &checkNodeN, mat_u64 &varNodeN, u64 first, u64 second)
{
    vec_u64 first_tmp = checkNodeN[first];
    vec_u64 second_tmp = checkNodeN[second];

    ldpc_code::zero_row(checkNodeN, varNodeN, first);
    ldpc_code::zero_row(checkNodeN, varNodeN, second);

    ldpc_code::add_rows(checkNodeN, varNodeN, first, second_tmp);
    ldpc_code::add_rows(checkNodeN, varNodeN, second, first_tmp);
}

void ldpc_code::swap_cols(mat_u64 &checkNodeN, mat_u64 &varNodeN, u64 first, u64 second)
{
    vec_u64 first_tmp = varNodeN[first];
    vec_u64 second_tmp = varNodeN[second];

    ldpc_code::zero_col(checkNodeN, varNodeN, first);
    ldpc_code::zero_col(checkNodeN, varNodeN, second);

    ldpc_code::add_cols(checkNodeN, varNodeN, first, second_tmp);
    ldpc_code::add_cols(checkNodeN, varNodeN, second, first_tmp);
}

void ldpc_code::add_rows(mat_u64 &checkNodeN, mat_u64 &varNodeN, u64 dest, const vec_u64 &src)
{
    vec_u64 new_row = checkNodeN[dest];
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

void ldpc_code::add_cols(mat_u64 &checkNodeN, mat_u64 &varNodeN, u64 dest, const vec_u64 &src)
{
    vec_u64 new_col = varNodeN[dest];
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

void ldpc_code::zero_row(mat_u64 &checkNodeN, mat_u64 &varNodeN, u64 m)
{
    for (auto vn : checkNodeN[m]) // from selected row, for each vn index, remove m from vn
    {
        varNodeN[vn].erase(std::remove(varNodeN[vn].begin(), varNodeN[vn].end(), m), varNodeN[vn].end());
    }
    checkNodeN[m] = vec_u64();
}

void ldpc_code::zero_col(mat_u64 &checkNodeN, mat_u64 &varNodeN, u64 n)
{
    for (auto cn : varNodeN[n]) // from selected col, for each cn index, remove n from cn
    {
        checkNodeN[cn].erase(std::remove(checkNodeN[cn].begin(), checkNodeN[cn].end(), n), checkNodeN[cn].end());
    }
    varNodeN[n] = vec_u64();
}


/**
 * @brief Prints parameters of LDPC code
 * 
 */
std::ostream &operator<<(std::ostream &os, const ldpc_code &code)
{
    std::cout << "nc : " << code.mN << "\n";
    std::cout << "mc : " << code.mM << "\n";
    std::cout << "kc : " << code.mK << "\n";
    std::cout << "nnz : " << code.mNNZ << "\n";
    std::cout << "nct :" << code.mNCT << "\n";
    std::cout << "mct : " << code.mMCT << "\n";
    std::cout << "kct : " << code.mKCT << "\n";
    std::cout << "max dc : " << code.mMaxDC << "\n";
    std::cout << "num puncture: " << code.mPuncture.size() << "\n";
    std::cout << "puncture: ";
    for (auto x : code.mPuncture)
    {
        std::cout << x << " ";
    }
    std::cout << "\nnum shorten: " << code.mShorten.size() << "\n";
    std::cout << "shorten: ";
    for (auto x : code.mShorten)
    {
        std::cout << x << " ";
    }
    return os;
}

} // namespace ldpc
