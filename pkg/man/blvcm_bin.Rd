\name{blvcm_bin}
\alias{blvcm_bin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bayesian latent variable collapsing model (BLVCM) for binary data with probit link
}
\description{
The function implements BLVCM for binary traits using a Gibbs sampler with probit link function. 
}
\usage{
blvcm_bin(pheno, geno, model = 3, iter = 30000, burnin = 500, var = -1, lambda = 0.2, 
cov = 0, init = c(0, 0))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pheno}{
An \eqn{N} x \eqn{3} phenotypic data matrix (trait, family number, zyg={1 for MZ, 2 for DZ}), where \eqn{N} is the number of subjects. The trait must be 0 or 1.
}
  \item{geno}{
An \eqn{N} x \eqn{K} genotypic data matrix, where \eqn{N} is the number of subjects and \eqn{K} is the number of rare variants. The value can be 0 or 1. A missing genotype is represented by -9, which will be imputated by BLVCM based on HWE.
}
  \item{model}{
Twin model: 3 for ACE model, 2 for AE model, 1 for independent subjects
}
  \item{iter}{
The number of MCMC iterations (must be positive). The default value is 30000.
}

 \item{burnin}{
The number of burn-ins (must be positive). The default value is 500.
}

  \item{var}{
The variance hyperparameters (must be positive) in the priors for \eqn{\beta} and \eqn{\gamma}. The default value is 1.
}
  \item{lambda}{
The threshold \eqn{\lambda} (must be positive) for hypothesis test. The default value is 0.2.
}
  \item{cov}{
A matrix of other covariates to be adjusted.
}
 \item{init}{
Initial values for \eqn{\beta} and \eqn{\gamma}. The default values are 0. The initial value for \eqn{\beta} must be non-negative.
}

}
\details{
The Gibbs sampler uses the variable augmentation method for probit link described in Albert, J. H., & Chib, S. (1993).  
Since the variance of a binary variable is determined by its mean compared to quantitative traits, \eqn{\theta}(s) are eliminated to avoid overfitting.

}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{BF_main }{The Bayes factor of the main effect}
  \item{BF_int }{The Bayes factor of the interaction effect}
  \item{post_odds_beta }{The posterior odds of \eqn{\beta}}
  \item{post_odds_gamma }{The posterior odds of \eqn{\gamma}}
  \item{com_a }{The inverse of the posterior mean of the precision for additive genetic component}
  \item{com_c }{The inverse of the posterior mean of the precision for shared environmental component}
  \item{mean_mu }{The posterior mean of the intercept \eqn{\mu}}
  \item{mean_beta }{The posterior mean of \eqn{\beta}}
  \item{mean_gamma }{The posterior mean of \eqn{\gamma}}
  \item{sd_mu }{The posterior standard deviation of the intercept \eqn{\mu}}
  \item{sd_beta }{The posterior standard deviation of \eqn{\beta}}
  \item{sd_gamma }{The posterior standard deviation of \eqn{\gamma}}
  \item{mean_rv }{The posterior mean of \eqn{\alpha}}
  \item{mean_cov }{The posterior mean of the effects of covariates}
%% ...
}
\references{
He, L., Sillanp‰‰, M. J., Ripatti, S., & Pitk‰niemi, J. (2014). Bayesian Latent Variable Collapsing Model for Detecting Rare Variant Interaction Effect in Twin Study. Genetic epidemiology, 38(4), 310-324.

Albert, J. H., & Chib, S. (1993). Bayesian analysis of binary and polychotomous response data. Journal of the American statistical Association, 88(422), 669-679.
}
\author{
Liang He
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
data(blvcm_bin_data)
blvcm_bin(blvcm_bin_data$pheno, blvcm_bin_data$geno[,1:3], iter=5000, burnin=500, model=2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
