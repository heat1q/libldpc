#pragma once

#include "functions.h"

namespace ldpc
{

    /**
 * @brief LDPC code class
 * 
 */
    class ldpc_code
    {
    public:
        ldpc_code(const std::string &pFileName);

        u64 calc_rank();

        friend std::ostream &operator<<(std::ostream &os, const ldpc_code &code);

        //getter functions
        u64 nc() const { return mN; };
        u64 kc() const { return mN-mM; };
        u64 mc() const { return mM; };
        u64 nnz() const { return mNNZ; };
        const mat_u64 &cn() const { return mCN; };
        const mat_u64 &vn() const { return mVN; };
        const vec_u64 &r() const { return mEdgeCN; };
        const vec_u64 &c() const { return mEdgeVN; };
        u64 nct() const { return mNCT; };
        u64 kct() const { return mNCT-mMCT; };
        u64 mct() const { return mMCT; };
        const vec_u64 &puncture() const { return mPuncture; };
        const vec_u64 &shorten() const { return mShorten; };
        u64 max_dc() const { return mMaxDC; };
        const mat_u64 &cn_index() const { return mCheckNodeN; };
        const mat_u64 &vn_index() const { return mVarNodeN; };

        const vec_u64 &bit_pos() const { return mBitPos; }

        // operations for gaussian elimination on sparse matrices over gf(2)
        static void swap_rows(mat_u64 &checkNodeN, mat_u64 &varNodeN, u64 i, u64 j);
        static void swap_cols(mat_u64 &checkNodeN, mat_u64 &varNodeN, u64 i, u64 j);
        static void add_rows(mat_u64 &checkNodeN, mat_u64 &varNodeN, u64 dest, const vec_u64 &src);
        static void add_cols(mat_u64 &checkNodeN, mat_u64 &varNodeN, u64 dest, const vec_u64 &src);
        static void zero_row(mat_u64 &checkNodeN, mat_u64 &varNodeN, u64 m);
        static void zero_col(mat_u64 &checkNodeN, mat_u64 &varNodeN, u64 n);

    private:
        u64 mN;
        u64 mM;
        u64 mNNZ;
        mat_u64 mCN;       /* denotes the check neighbors, i.e. connected VN, for each check node as index in c/r; dimensions cn[mc][cw[i]] */
        mat_u64 mVN;       /* denotes the var neighbors, i.e., connected CN, for each variable node as index in c/r; dimensions vn[nc][vw[i]] */
        vec_u64 mEdgeCN;   /* edge oriented non zero row indices; length nnz */
        vec_u64 mEdgeVN;   /* edge oriented non zero column indices; length nnz */
        vec_u64 mPuncture; /* array pf punctured bit indices */
        vec_u64 mShorten;  /* array of shortened bit indices */
        u64 mNCT;          /* number of transmitted code bits */
        u64 mMCT;          /* number of transmitted parity check bits */
        u64 mMaxDC;

        double mRate;

        // position of transmitted bits, i.e. puncture/shorten exluded
        vec_u64 mBitPos;

        mat_u64 mCheckNodeN; // direct check node neighbourhood with node oriented index
        mat_u64 mVarNodeN;   // direct variable node neighbourhood with node oriented index
    };

} // namespace ldpc
