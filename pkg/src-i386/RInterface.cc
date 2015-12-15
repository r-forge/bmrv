#include "paper_3.h"
#include "paper_2.h"
#include <R.h> // R functions
// #include "Rmath.h" // Rmath

//
void CppGibbs_hbmr(double *result, int * numRows, int * numCols, int * numCols2, int * fam, double *mat1, double *mat2, double * mat3, double * mat4, double * mat_kin, double * arg_m, double * arg_i, int * arg_n, int * arg_bu, double * arg_t, int * gamma_est, double * arg_a, double * arg_b)
{
	int rvinfo = 0;
//
	if((*gamma_est)!=0)
	{
		rvinfo = 1;
	}

	gibbssampler2(result, numRows, numCols, numCols2, fam, mat1, mat2, mat3, mat4, mat_kin, arg_m, arg_i, arg_n, arg_bu, arg_t, arg_a, arg_b, rvinfo);

}
//-----------------------------------------------------------------------------

void CppGibbs_hbmr_bin(double *result, int * numRows, int * numCols, int * numCols2, int * fam, int *mat1, double *mat2, double * mat3, double * mat4, double * mat_kin, double * arg_m, double * arg_i, int * arg_n, int * arg_bu, double * arg_t, int * gamma_est, double * arg_a, double * arg_b)
{

//
	gibbssampler2_bin(result, numRows, numCols, numCols2, fam, mat1, mat2, mat3, mat4, mat_kin, arg_m, arg_i, arg_n, arg_bu, arg_t, arg_a, arg_b);

}

void CppGibbs_hbmr_ord(double *result, int * numRows, int * numCols, int * numCols2, int * fam, int *mat1, double *mat2, double * mat3, double * mat4, double * mat_kin, double * arg_m, double * arg_i, int * arg_n, int * arg_bu, double * arg_t, int * gamma_est, double * arg_a, double * arg_b)
{

//
	gibbssampler2_ord(result, numRows, numCols, numCols2, fam, mat1, mat2, mat3, mat4, mat_kin, arg_m, arg_i, arg_n, arg_bu, arg_t, arg_a, arg_b);

}

void CppGibbs_blvcm(double *result,int * numRows, int * numCols, int * numCols2, double *mat1, int *mat2, double * mat3, double * arg_v, double * arg_i, int * arg_n, int * arg_b, int * arg_t, double * arg_i_b, double * arg_i_g)
{

gibbssampler(result,numRows, numCols, numCols2,mat1, mat2, mat3, arg_v, arg_i, arg_n, arg_b, arg_t, arg_i_b, arg_i_g);

}

void CppGibbs_blvcm_bin(double *result,int * numRows, int * numCols, int * numCols2,int *mat1, int *mat2, double * mat3, double * arg_v, double * arg_i, int * arg_n, int * arg_b, int * arg_t, double * arg_i_b, double * arg_i_g)
{

gibbssampler_bin(result,numRows, numCols, numCols2,mat1, mat2, mat3, arg_v, arg_i, arg_n, arg_b, arg_t, arg_i_b, arg_i_g);

}

//
extern "C" {
void CWrapper_hbmr(double *result, int * numRows, int * numCols, int * numCols2, int * fam, double *mat1, double *mat2, double * mat3, double * mat4, double * mat_kin, double * arg_m, double * arg_i, int * arg_n, int * arg_bu, double * arg_t, int * gamma_est, double * arg_a, double * arg_b)
{
//- Invoke second function which internally can do C++ things.
//
CppGibbs_hbmr(result, numRows, numCols, numCols2, fam, mat1, mat2, mat3, mat4, mat_kin, arg_m, arg_i, arg_n, arg_bu, arg_t, gamma_est, arg_a, arg_b);
}

// HBMR for binary traits
void CWrapper_hbmr_bin(double *result, int * numRows, int * numCols, int * numCols2, int * fam, int *mat1, double *mat2, double * mat3, double * mat4, double * mat_kin, double * arg_m, double * arg_i, int * arg_n, int * arg_bu, double * arg_t, int * gamma_est, double * arg_a, double * arg_b)
{
//- Invoke second function which internally can do C++ things.
//
CppGibbs_hbmr_bin(result, numRows, numCols, numCols2, fam, mat1, mat2, mat3, mat4, mat_kin, arg_m, arg_i, arg_n, arg_bu, arg_t, gamma_est, arg_a, arg_b);
}


// HBMR for ordinal traits
void CWrapper_hbmr_ord(double *result, int * numRows, int * numCols, int * numCols2, int * fam, int *mat1, double *mat2, double * mat3, double * mat4, double * mat_kin, double * arg_m, double * arg_i, int * arg_n, int * arg_bu, double * arg_t, int * gamma_est, double * arg_a, double * arg_b)
{
//- Invoke second function which internally can do C++ things.
//
CppGibbs_hbmr_ord(result, numRows, numCols, numCols2, fam, mat1, mat2, mat3, mat4, mat_kin, arg_m, arg_i, arg_n, arg_bu, arg_t, gamma_est, arg_a, arg_b);
}


void CWrapper_blvcm(double *result,int * numRows, int * numCols, int * numCols2,double *mat1, int *mat2, double * mat3, double * arg_v, double * arg_i, int * arg_n, int * arg_b, int * arg_t, double * arg_i_b, double * arg_i_g)
{
//- Invoke second function which internally can do C++ things.
//
CppGibbs_blvcm(result,numRows, numCols, numCols2,mat1, mat2, mat3, arg_v, arg_i, arg_n, arg_b, arg_t, arg_i_b, arg_i_g);
}


// BLVCM for binary traits
void CWrapper_blvcm_bin(double *result,int * numRows, int * numCols, int * numCols2,int *mat1, int *mat2, double * mat3, double * arg_v, double * arg_i, int * arg_n, int * arg_b, int * arg_t, double * arg_i_b, double * arg_i_g)
{
//- Invoke second function which internally can do C++ things.
//
CppGibbs_blvcm_bin(result,numRows, numCols, numCols2,mat1, mat2, mat3, arg_v, arg_i, arg_n, arg_b, arg_t, arg_i_b, arg_i_g);
}
}