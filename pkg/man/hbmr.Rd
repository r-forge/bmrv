\name{hbmr}
\alias{hbmr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hierarchical Bayesian multiple regression model incorporating genotype uncertainty (HBMR)
}
\description{
The function implements HBMR using Gibbs sampling method for quantitative traits.
}
\usage{
hbmr(pheno, geno, qi = matrix(), fam = 0, kin = matrix(), iter = 10000, burnin = 500, 
gq = 20, imp = 0.1, cov = matrix(), maf = c(), rvinfo = FALSE, pa = 1.3, pb = 0.04)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pheno}{
Phenotypic vector (\eqn{N} x \eqn{1}). For better inference, it is recommanded that phenotype should be standardized.}
  \item{geno}{
\eqn{N} x \eqn{K} Genotypic data matrix, where \eqn{N} is the number of subjects and \eqn{K} is the number of rare variants.
	Genotypic value is only for dominant coding, i.e. 0 or 1. Plug in 0 for imputed genotypes.
}
  \item{qi}{ 
An optional \eqn{N} x \eqn{K} Genotypic quality matrix, where \eqn{N} is the number of subjects and \eqn{K} is the number of rare variants.
	If the genotype is sequenced, this must be an integer >=1 and is its \eqn{GQ} score in VCF file. If the genotype is imputed, this must be a value <1, and is its expected genotypic value based on the dominant coding. 
}
  \item{fam}{ 
  fam=1 for family samples. In this case, a relatedness matrix should be given. See kin.
}
  \item{kin}{ 
  In the case of fam=1, kin is an \eqn{N} x \eqn{N} relatedness matrix. 
}
  \item{iter}{
The number of MCMC iterations. The default value is 10000.
}

 \item{burnin}{
The number of burn-ins. The default value is 500.
}

  \item{gq}{
A cutoff for \eqn{GQ} score (\eqn{\lambda_Q}). It should be an positive integer. If not specified, default value is 20. See the reference for more details.
}
  \item{imp}{
A cutoff for imputed genotype (\eqn{\lambda_I}). It should be a real number in (0,1). If not specified, default value is 0.1. See the reference for more details.
}
  \item{cov}{
An optional \eqn{N} x \eqn{M} covariate data matrix, where N is the number of subjects and M is the number of covariates.
}
  \item{maf}{
An optional minor allele frequency information vector (\eqn{K} x 1). If not specified, MAF will be estimated based on the genotype data.
}
  \item{rvinfo}{
TRUE or FALSE. Default is FALSE. Indicator of showing estimatd RV effect size and standard deviation.
}
  \item{pa}{
The positive hyper-parameter \eqn{a} in the gamma distribution of Bayesian shrinkage prior. The default value is 1.3.
}
  \item{pb}{
The positive hyper-parameter \eqn{b} in the gamma distribution of Bayesian shrinkage prior. The default value is 0.04.
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{BF }{The Bayes factor of \eqn{\delta=1} vs. \eqn{\delta=0}}
  \item{BF_RB }{The BF estimated by using Rao-Blackwellization theorem}
  \item{p_upper }{For a BF larger than 2, we calculate p_upper that is the upper bound of the p value corresponding to the BF based on the connection \eqn{BF<(-1)/(e*p*log(p))}. The exact p value, which is smaller than p_upper, can be obtained through permutations.}
  \item{mean }{The mean of the posterior of \eqn{\beta_0}}
  \item{var }{The inverse of the mean of posterior of precision 1/\eqn{\sigma}}
  \item{est_geno }{The number of genotypes whose uncertainty are considered in estimation}
  \item{var_ran }{The estimated variance of the random effect for family design}
  \item{rv_mean_es }{The means of the posterior of \eqn{\gamma} for the \eqn{K} RVs}
  \item{rv_sd_es }{The standard deviations of the posterior of \eqn{\gamma} for the \eqn{K} RVs}
  \item{mean_cov }{The means of the posterior of for the \eqn{M} covariates}
%% ...
}
\references{
He, L., Pitkäniemi, J., Sarin, A. P., Salomaa, V., Sillanpää, M. J., & Ripatti, S. (2015). Hierarchical Bayesian Model for Rare Variant Association Analysis Integrating Genotype Uncertainty in Human Sequence Data. Genetic epidemiology, 39(2), 89-100.
}
\author{
Liang He
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
data(hbmr_data)
hbmr(hbmr_data$pheno_data, hbmr_data$geno_data[,1:3], hbmr_data$qual_data[,1:3], 
iter=10000, burnin=1000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
