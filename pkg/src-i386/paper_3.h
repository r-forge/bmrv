#include "stdlib.h"
#include "math.h"

#include <fstream>
#include <iostream>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>

#define THRES_QI 20
#define MAX_RV 1000
// #define AFO 99
#define MAX_COV 1000
// #define INTERC -2

void gibbssampler2(double *result, int * numRows, int * numCols, int * numCols2, int * fam, double *mat1, double *mat2, double * mat3, double * mat4, double * mat_kin, double * arg_m, double * arg_i, int * arg_n, int * arg_bu, double * arg_t, double * arg_a, double * arg_b, int rvinfo);

void gibbssampler2_bin(double *result, int * numRows, int * numCols, int * numCols2, int * fam, int *mat1, double *mat2, double * mat3, double * mat4, double * mat_kin, double * arg_m, double * arg_i, int * arg_n, int * arg_bu, double * arg_t, double * arg_a, double * arg_b);

void gibbssampler2_ord(double *result, int * numRows, int * numCols, int * numCols2, int * fam, int *mat1, double *mat2, double * mat3, double * mat4, double * mat_kin, double * arg_m, double * arg_i, int * arg_n, int * arg_bu, double * arg_t, double * arg_a, double * arg_b);
