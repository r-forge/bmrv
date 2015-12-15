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
#include <boost/random/exponential_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>

#define MAX_RV 1000
// #define AFO 99
#define MAX_COV 1000
// #define INTERC -2


void gibbssampler(double *result,int * numRows, int * numCols, int * numCols2,double *mat1, int *mat2, double * mat3, double * arg_v, double * arg_i, int * arg_n, int * arg_b, int * arg_t, double * arg_i_b, double * arg_i_g);

void gibbssampler_bin(double *result,int * numRows, int * numCols, int * numCols2, int *mat1, int *mat2, double * mat3, double * arg_v, double * arg_i, int * arg_n, int * arg_b, int * arg_t, double * arg_i_b, double * arg_i_g);
