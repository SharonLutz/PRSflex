\name{PRSflex}
\alias{PRSflex}
\title{PRSflex}
\description{
This needs filled in.
}
\usage{
PRSflex(nSim = 100, nTrain = 2500, nTest = 2500, varY = 1, sigl = 5e-08, betaS = seq(0, 0.1, by = 0.01), nAdd = 25, nDom = 25, nCod = 25, nRec = 25, MAFadd = 0.05, MAFdom = 0.05, MAFcod = 0.05, MAFrec = 0.05, plot.pdf = TRUE, SEED = 1)
}
\arguments{
  \item{nSim}{
  is the number of simulations run for the function.
}
  \item{nTrain}{
  is the sample size of the training dataset.
}
  \item{nTest}{
  is the sample size of the test dataset.
}
  \item{varY}{
  is the variance for outcome Y.
}
  \item{sigl}{
  is the significance level.
}
  \item{betaS}{
  is a sequence of the effect sizes.
}
  \item{nAdd}{
  is the number of additive SNPs.
}
  \item{nDom}{
  is the number of dominant SNPs.
}
  \item{nCod}{
  is the number of codominant SNPs.
}
  \item{nRec}{
  is the number of recessive SNPs.
}
  \item{MAFadd}{
  is the minor allele frequency for additive SNPs.
}
  \item{MAFdom}{
  is the minor allele frequency for dominant SNPs.
}
  \item{MAFcod}{
  is the minor allele frequency for codominant SNPs.
}
  \item{MAFrec}{
  is the minor allele frequency for recessive SNPs.
}
  \item{plot.pdf}{
  saves three PDF plots of results to working directory if set to TRUE. If FALSE, plots are not saved. 
}
  \item{SEED}{
  sets the seed.
}
}

\value{
The function outputs a matrix with columns for average mean squared error, average adjusted R^2, and average proportion of signicant p-values for each method and each effect size.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Sharon Lutz, Kirsten Voorhies
}
\section{Warning}{
  library(caret), library(stats), and library(HDeconometrics) are needed to run this function.
}


\examples{


PRSflex(nSim = 50, nTrain = 1000, nTest = 1000, varY = 1, 
    sigl = 5e-08, betaS = seq(0, 0.08, by = 0.01), nAdd = 25, 
    nDom = 25, nCod = 25, nRec = 25, MAFadd = 0.05, MAFdom = 0.05, 
    MAFcod = 0.05, MAFrec = 0.05, plot.pdf = FALSE, SEED = 1) 

}

\keyword{ PRSflex }

