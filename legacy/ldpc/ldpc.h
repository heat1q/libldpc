#pragma once

#include <iostream>
#include <vector>
#include <memory>
#include <cmath>
#include <cstring>
#include <chrono>
#include <fstream>

namespace ldpc
{
using bits_t = unsigned char;
using labels_t = unsigned short;
using symbols_t = unsigned;

using vec_bits_t = std::vector<bits_t>;
using vec_labels_t = std::vector<unsigned short>;
using vec_symbols_t = std::vector<unsigned>;
using vec_size_t = std::vector<std::size_t>;
using vec_double_t = std::vector<double>;

using mat_bits_t = std::vector<std::vector<bits_t>>;
using mat_size_t = std::vector<std::vector<std::size_t>>;
using mat_double_t = std::vector<std::vector<double>>;

class ldpc_code
{
public:
	ldpc_code(const char *pFileName, const char *pClFile);
	void print();

	//getter functions
	std::size_t nc() const { return mN; };
	std::size_t kc() const { return mK; };
	std::size_t mc() const { return mM; };
	std::size_t nnz() const { return mNNZ; };
	const vec_size_t &cw() const { return mCW; };
	const vec_size_t &vw() const { return mVW; };
	const mat_size_t &cn() const { return mCN; };
	const mat_size_t &vn() const { return mVN; };
	const vec_size_t &r() const { return mR; };
	const vec_size_t &c() const { return mC; };
	std::size_t nct() const { return mNCT; };
	std::size_t kct() const { return mKCT; };
	std::size_t mct() const { return mMCT; };
	const vec_size_t &puncture() const { return mPuncture; };
	std::size_t num_puncture() const { return mNumPuncture; };
	const vec_size_t &shorten() const { return mShorten; };
	std::size_t num_shorten() const { return mNumShorten; };
	std::size_t max_dc() const { return mMaxDC; };
	std::size_t nl() const { return mNL; };
	const mat_size_t &layers() const { return mLayers; };

private:
	std::size_t mN;
	std::size_t mK;
	std::size_t mM;
	std::size_t mNNZ;
	vec_size_t mCW;				 /* denotes the check weight of each check node, i.e., # of connected VN; dimensions cw[mc] */
	vec_size_t mVW;				 /* denotes the variable weight, i.e., # of connected CN; dimensions vw[nc] */
	mat_size_t mCN;				 /* denotes the check neighbors, i.e. connected VN, for each check node as index in c/r; dimensions cn[mc][cw[i]] */
	mat_size_t mVN;				 /* denotes the var neighbors, i.e., connected CN, for each variable node as index in c/r; dimensions vn[nc][vw[i]] */
	vec_size_t mR;				 /* non zero row indices; length nnz */
	vec_size_t mC;				 /* non zero check indices; length nnz */
	vec_size_t mPuncture;		 /* array pf punctured bit indices */
	std::size_t mNumPuncture;	/* number of punctured bits */
	std::size_t mNumPunctureSys; /* number of punctured bits in systematic part */
	std::size_t mNumPuncturePar; /* number of punctured bits in parity part */
	vec_size_t mShorten;		 /* array of shortened bit indices */
	std::size_t mNumShorten;	 /* number of shortened bits */
	std::size_t mNCT;			 /* number of transmitted code bits */
	std::size_t mKCT;			 /* number of transmitted information bits */
	std::size_t mMCT;			 /* number of transmitted parity check bits */
	std::size_t mMaxDC;
	std::size_t mNL;
	mat_size_t mLayers;
};

class ldpc_decoder
{
public:
	ldpc_decoder(ldpc_code *pCode, std::size_t pI, bool pEarlyTerm);

	std::size_t decode_layered();
	bool is_codeword_legacy();

	std::size_t max_iter() const { return mMaxIter; }
	bool early_termination() const { return mEarlyTerm; }

	ldpc_code* mLdpcCode;

	vec_double_t mLv2c;
	vec_double_t mLc2v;
	vec_double_t mExMsgCN;

	vec_double_t mLLRIn;
	vec_double_t mLLROut;

	vec_bits_t mSynd;
	vec_bits_t mCO;

	std::size_t mIter;

	bool mIsCW;

	std::size_t mMaxIter;
	bool mEarlyTerm;
};

void dec2bin(std::size_t val, uint8_t m);
double jacobian(double L1, double L2);
double jacobian_lin_approx(double L);
int sign(double a);

} // namespace ldpc