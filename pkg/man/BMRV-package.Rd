\name{BMRV-package}
\alias{BMRV-package}
\alias{BMRV}
\docType{package}
\title{
Bayesian models for rare variant association analysis
}
\description{
This package provides two Bayesian models for detecting the association between a set of rare genetic variants and a trait that can be continuous, ordinal or binary. 
BLVCM detects interaction effect and is dedicated to twin design while it can also be applied to independent samples. 
HBMR incorporates genotype uncertainty information and can be applied to either independent or family samples. Furhtermore, it deals with continuous, binary and ordinal traits.
}
\details{
\tabular{ll}{
Package: \tab BMRV\cr
Type: \tab Package\cr
Version: \tab 1.32\cr
Date: \tab 2016-10-29\cr
License: \tab GPL (>=2)\cr
}
blvcm
hbmr
}
\author{
Liang He

Maintainer: Liang He <liang.he@duke.edu>
}
\references{
He, L., Sillanp‰‰, M. J., Ripatti, S., & Pitk‰niemi, J. (2014). Bayesian Latent Variable Collapsing Model for Detecting Rare Variant Interaction Effect in Twin Study. Genetic epidemiology, 38(4), 310-324.

He, L., Pitk‰niemi, J., Sarin, A. P., Salomaa, V., Sillanp‰‰, M. J., & Ripatti, S. (2015). Hierarchical Bayesian Model for Rare Variant Association Analysis Integrating Genotype Uncertainty in Human Sequence Data. Genetic epidemiology, 39(2), 89-100.
}

\keyword{ package }

\examples{
data(blvcm_data)
temp<- blvcm(blvcm_data$pheno_data, blvcm_data$geno_data[,1:3], iter=10000, model = 3)
}
