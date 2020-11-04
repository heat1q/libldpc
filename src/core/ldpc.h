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
        /**
         * @brief Construct a new ldpc code object.
         * 
         * @param pcFileName parity-check matrix file
         */
        ldpc_code(const std::string &pcFileName);

        /**
         * @brief Construct a new ldpc code object
         * 
         * @param pcFileName parity-check matrix file
         * @param genFileName generator matrix file
         */
        ldpc_code(const std::string &pcFileName, const std::string &genFileName);

        /**
         * @brief Read the parity-check matrix from file.
         * 
         * @param pcFileName Filename
         */
        void read_H(const std::string &pcFileName);

        /**
         * @brief Read the generator matrix from file.
         * 
         * @param genFileName Filename
         */
        void read_G(const std::string &genFileName);

        /**
         * @brief Calculates the rank of H.
         * 
         * @return u64 rank
         */
        u64 calc_rank();

        friend std::ostream &operator<<(std::ostream &os, const ldpc_code &code);

        // Number of columns (variable nodes)
        int nc() const { return mH.num_cols(); };
        // Number of information bits
        int kc() const { return mH.num_cols() - mH.num_rows(); };
        // Number of parity-checks (check nodes)
        int mc() const { return mH.num_rows(); };
        // Number of non-zero entries (number of edges)
        int nnz() const { return mH.nz_entry().size(); };
        // Number of transmitted columns (variable nodes)
        int nct() const { return nc() - mPuncture.size() - mShorten.size(); };
        // Number of transmitted information bits
        int kct() const { return nct() - mct(); };
        // Number of transmitted parity-checks (check nodes)
        int mct() const { return mc() - mPuncture.size(); };
        // Array of puncture indices
        const vec_u64 &puncture() const { return mPuncture; };
        // Array of shorten indices
        const vec_u64 &shorten() const { return mShorten; };
        // Maximum check node degree
        u64 max_dc() const { return mMaxDC; };
        const vec_u64 &bit_pos() const { return mBitPos; }
        // Parity-check matrix
        const sparse_csr<bits_t> &H() const { return mH; }
        // Generator matrix
        const sparse_csr<bits_t> &G() const { return mG; }

        // operations for gaussian elimination on sparse matrices over gf(2)
        static void swap_rows(mat_u64 &checkNodeN, mat_u64 &varNodeN, u64 i, u64 j);
        static void swap_cols(mat_u64 &checkNodeN, mat_u64 &varNodeN, u64 i, u64 j);
        static void add_rows(mat_u64 &checkNodeN, mat_u64 &varNodeN, u64 dest, const vec_u64 &src);
        static void add_cols(mat_u64 &checkNodeN, mat_u64 &varNodeN, u64 dest, const vec_u64 &src);
        static void zero_row(mat_u64 &checkNodeN, mat_u64 &varNodeN, u64 m);
        static void zero_col(mat_u64 &checkNodeN, mat_u64 &varNodeN, u64 n);

        
    private:
        vec_u64 mPuncture; /* array pf punctured bit indices */
        vec_u64 mShorten;  /* array of shortened bit indices */
        u64 mMaxDC;

        // position of transmitted bits, i.e. puncture/shorten exluded
        vec_u64 mBitPos;

        sparse_csr<bits_t> mH; // Parity-Check Matrix
        sparse_csr<bits_t> mG; // Generator Matrix
    };

} // namespace ldpc
