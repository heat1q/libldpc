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
        const vec_int &puncture() const { return mPuncture; };
        // Array of shorten indices
        const vec_int &shorten() const { return mShorten; };
        // Maximum check node degree
        int max_dc() const { return mMaxDC; };
        // Index position of transmitted bits 
        const vec_int &bit_pos() const { return mBitPos; }
        // Parity-check matrix
        const sparse_csr<bits_t> &H() const { return mH; }
        // Generator matrix
        const sparse_csr<bits_t> &G() const { return mG; }
        
    private:
        vec_int mPuncture; /* array pf punctured bit indices */
        vec_int mShorten;  /* array of shortened bit indices */
        int mMaxDC;

        // position of transmitted bits, i.e. puncture/shorten exluded
        vec_int mBitPos;

        sparse_csr<bits_t> mH; // Parity-Check Matrix
        sparse_csr<bits_t> mG; // Generator Matrix
    };

} // namespace ldpc
