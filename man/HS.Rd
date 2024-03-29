% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HS.R
\name{HS}
\alias{HS}
\title{Hierarchical Spike}
\usage{
HS(
  Betah,
  Sigmah,
  kappa0 = kappa0,
  kappastar0 = kappastar0,
  sigma20 = sigma20,
  s20 = s20,
  m,
  K,
  niter = 1000,
  burnin = 500,
  nthin = 2,
  nchains = 2,
  a1 = 0.1,
  a2 = 0.1,
  d1 = 0.1,
  d2 = 0.1,
  c1 = 1,
  c2 = 1,
  e2 = 1,
  snpnames,
  genename
)
}
\arguments{
\item{Betah}{A list containing m-dimensional vectors of the regression coefficients for K studies.}

\item{Sigmah}{A list containing the positive definite covariance matrices (m*m-dimensional) which is the estimated covariance matrices of K studies.}

\item{kappa0}{Initial value for kappa (its dimension is equal to nchains).}

\item{kappastar0}{Initial value for kappastar (its dimension is equal to nchains).}

\item{sigma20}{Initial value for sigma2 (its dimension is equal to nchains).}

\item{s20}{Initial value for s2 (its dimension is equal to nchains).}

\item{m}{Number of variables in the group.}

\item{K}{Number of traits.}

\item{niter}{Number of iterations for the Gibbs sampler.}

\item{burnin}{Number of burn-in iterations.}

\item{nthin}{The lag of the iterations used for the posterior analysis is defined (or thinning rate).}

\item{nchains}{Number of Markov chains, when nchains>1, the function calculates the Gelman-Rubin convergence statistic, as modified by Brooks and Gelman (1998).}

\item{a1, a2}{Hyperparameters of kappa. Default is a1=0.1 and a2=0.1.}

\item{d1, d2}{Hyperparameters of sigma2. Default is d1=0.1 and d2=0.1.}

\item{c1, c2}{Hyperparameters of kappastar. Default is c1=0.1 and c2=0.1.}

\item{e2}{Initial value for doing Monte Carlo EM algorithm to estimate hyperparameter of s2.}

\item{snpnames}{Names of variables for the group.}

\item{genename}{Name of group.}
}
\value{
\itemize{
\item mcmcchain: The list of simulation output for all parameters.
\item Summary: Summary statistics for regression coefficients in each study.
\item Indicator 1: A table containing m rows of binary indicators for each study, the number of studies with nonzero signal and having pleiotropic effect by credible interval (CI). The first K columns show nonzero signals, K+1 th column includes the number of studies with nonzero signal and the last column shows an indicator for having pleiotropic effect of each SNP.
\item Indicator 2: A table containing m rows of binary indicators for each study, the number of studies with nonzero signal and having pleiotropic effect by median. The first K columns show nonzero signals, K+1 th column includes the number of studies with nonzero signal and the last column shows an indicator for having pleiotropic effect of each SNP.
}
}
\description{
Run a Gibbs sampler for a multivariate Bayesian sparse group selection model with hierarchical spike prior for detecting pleiotropic effects on the traits. This function is designed for summary statistics containing estimated regression coefficients and their estimated covariance matrices.
}
\details{
For considering the HS prior, a reparameterization of betah_k is considered as betah_k = V_k^0.5 b_k,  V_k^0.5 = diag(tau_k1,...,tau_km)$. Therfore, we have the following hirarchical model:

b_k ~ (1 - kappa) delta_0(b_k) + kappa N_m(0,sigma2 I_m ),

tau_kj ~ (1 - kappa^*) delta_0(tau_kj) + kappa TN(0,s2),

kappa ~ Beta(a_1,a_2),

kappa^* ~ Beta(c_1,c_2),

sigma2 ~ inverseGamma (d_1,d_2).

s2 ~ inverseGamma (e_1,e_2).

where delta_0(betah_k) denotes a point mass at 0, such that delta_0(betah_k)=1 if beta_k=0 and  delta_0(betah_k)=0 if  at least one of the $m$ components of beta_k is non-zero and TN(0,s2)  denotes a univariate truncated normal distribution at zero with mean 0 and variance s2.
}
\examples{
############################# PARP2_summary ###############################################
data(PARP2_summary)
Breast <- PARP2_summary$Breast
Thyroid <- PARP2_summary$Thyroid
Betah <- list(Breast$beta, Thyroid$beta)
Sigmah <- list(diag(Breast$se), diag(Thyroid$se))
genename <- "PARP2"
snpnames <- Breast$snp
K <- 2
m <- 6

RES <- HS(Betah, Sigmah,
  kappa0 = 0.5, kappastar0 = 0.5, sigma20 = 1, s20 = 1,
  m = m, K = K, niter = 1000, burnin = 500, nthin = 1, nchains = 1,
  a1 = 0.1, a2 = 0.1, d1 = 0.1, d2 = 0.1, c1 = 1, c2 = 1, e2 = 1, snpnames, genename
)
\dontrun{
  print(RES)
  ############################# Gene DNAJC1 ###############################################
  data(DNAJC1)
  Breast <- DNAJC1$Breast
  Thyroid <- DNAJC1$Thyroid
  genename <- "DNAJC1"
  snpnames <- Breast$snp
  Betah <- list(Breast$beta, Thyroid$beta)
  Sigmah <- list(diag(Breast$se^2), diag(Thyroid$se^2))
  K <- 2
  m <- 14

  RES <- HS(Betah, Sigmah,
    kappa0 = 0.5, kappastar0 = 0.5, sigma20 = 1, s20 = 1,
    m = m, K = K, niter = 2000, burnin = 1000, nthin = 1, nchains = 1,
    a1 = 0.1, a2 = 0.1, d1 = 0.1, d2 = 0.1, c1 = 1, c2 = 1, e2 = 1, snpnames, genename
  )
  print(RES)

  RES1 <- HS(Betah, Sigmah,
    kappa0 = c(0.5, 0.3), kappastar0 = c(0.5, 0.3), sigma20 = c(2, 1), s20 = c(1, 2),
    m = m, K = K, niter = 2000, burnin = 1000, nthin = 2, nchains = 2,
    a1 = 0.1, a2 = 0.1, d1 = 0.1, d2 = 0.1, c1 = 1, c2 = 1, e2 = 1, snpnames, genename
  )
  print(RES1)
  ################### Simulated summary level data with K=5 ###############################
  data(Simulated_summary)
  genename <- Simulated_summary$genename
  snpnames <- Simulated_summary$snpnames
  Betah <- Simulated_summary$simBeta
  Sigmah <- Simulated_summary$simSIGMA
  K <- 5
  m <- 10

  RES <- HS(Betah, Sigmah,
    kappa0 = 0.5, kappastar0 = 0.5, sigma20 = 1, s20 = 1,
    m = m, K = K, niter = 2000, burnin = 1000, nthin = 2, nchains = 1,
    a1 = 0.1, a2 = 0.1, d1 = 0.1, d2 = 0.1, c1 = 1, c2 = 1, e2 = 1, snpnames, genename
  )
  print(RES)

  RES1 <- HS(Betah, Sigmah,
    kappa0 = c(0.5, 0.3), kappastar0 = c(0.5, 0.3), sigma20 = c(2, 1), s20 = c(1, 2),
    m = m, K = K, niter = 2000, burnin = 1000, nthin = 2, nchains = 2,
    a1 = 0.1, a2 = 0.1, d1 = 0.1, d2 = 0.1, c1 = 1, c2 = 1, e2 = 1, snpnames, genename
  )
  print(RES1)
  ################################### Gene PARP2 ##########################################
  library(BhGLM)
  data(PARP2)
  Breast <- PARP2$Breast
  Thyroid <- PARP2$Thyroid
  genename <- "PARP2"
  snpnames <- c("rs3093872", "rs3093921", "rs1713411", "rs3093926", "rs3093930", "rs878156")


  Fit1 <- BhGLM::bglm(y1 ~ ., family = binomial(link = "logit"), data = Breast)
  Betah1 <- Fit1$coefficients[-1]
  Sigmah1 <- cov(coef(arm::sim(Fit1)))[-1, -1]


  Fit2 <- BhGLM::bglm(y2 ~ ., family = binomial(link = "logit"), data = Thyroid)
  Betah2 <- Fit2$coefficients[-1]
  Sigmah2 <- cov(coef(arm::sim(Fit2)))[-1, -1]

  Betah <- list(Betah1, Betah2)
  Sigmah <- list(Sigmah1, Sigmah2)
  K <- 2
  m <- 6


  RES <- HS(Betah, Sigmah,
    kappa0 = 0.5, kappastar0 = 0.5, sigma20 = 1, s20 = 1,
    m = m, K = K, niter = 2000, burnin = 1000, nthin = 2, nchains = 1,
    a1 = 0.1, a2 = 0.1, d1 = 0.1, d2 = 0.1, c1 = 1, c2 = 1, e2 = 1, snpnames, genename
  )
  print(RES)

  RES1 <- HS(Betah, Sigmah,
    kappa0 = c(0.5, 0.3), kappastar0 = c(0.5, 0.3), sigma20 = c(2, 1), s20 = c(1, 2),
    m = m, K = K, niter = 2000, burnin = 1000, nthin = 2, nchains = 2,
    a1 = 0.1, a2 = 0.1, d1 = 0.1, d2 = 0.1, c1 = 1, c2 = 1, e2 = 1, snpnames, genename
  )
  print(RES1)
  ########### Simulated individual level data with K=3 and continuous phynotype ###########
  library(BhGLM)
  data(Simulated_individual)
  Study1 <- Simulated_individual$Study1
  Study2 <- Simulated_individual$Study2
  Study3 <- Simulated_individual$Study3
  K <- 3
  m <- 30
  genename <- "Simulated"
  snpnames <- sprintf("SNP\%s", seq(1:m))


  Fit1 <- BhGLM::bglm(Y1 ~ ., family = gaussian, data = data.frame(Study1))
  Betah1 <- Fit1$coefficients[-1]
  Sigmah1 <- cov(coef(arm::sim(Fit1)))[-1, -1]


  Fit2 <- BhGLM::bglm(Y2 ~ ., family = gaussian, data = data.frame(Study2))
  Betah2 <- Fit2$coefficients[-1]
  Sigmah2 <- cov(coef(arm::sim(Fit2)))[-1, -1]


  Fit3 <- BhGLM::bglm(Y3 ~ ., family = gaussian, data = data.frame(Study3))
  Betah3 <- Fit3$coefficients[-1]
  Sigmah3 <- cov(coef(arm::sim(Fit3)))[-1, -1]

  Betah <- list(Betah1, Betah2, Betah3)
  Sigmah <- list(Sigmah1, Sigmah2, Sigmah3)


  RES <- HS(Betah, Sigmah,
    kappa0 = 0.5, kappastar0 = 0.5, sigma20 = 1, s20 = 1,
    m = m, K = K, niter = 2000, burnin = 1000, nthin = 2, nchains = 1,
    a1 = 0.1, a2 = 0.1, d1 = 0.1, d2 = 0.1, c1 = 1, c2 = 1, e2 = 1, snpnames, genename
  )
  print(RES)

  RES1 <- HS(Betah, Sigmah,
    kappa0 = c(0.5, 0.3), kappastar0 = c(0.5, 0.3), sigma20 = c(2, 1), s20 = c(1, 2),
    m = m, K = K, niter = 2000, burnin = 1000, nthin = 2, nchains = 2,
    a1 = 0.1, a2 = 0.1, d1 = 0.1, d2 = 0.1, c1 = 1, c2 = 1, e2 = 1, snpnames, genename
  )
  print(RES1)

  ########### Simulated individual level data with K=2 and gene expression data ###########
  library(BhGLM)
  data(Simulated_individual_survival)
  Study1 <- Simulated_individual_survival$Study1
  Study2 <- Simulated_individual_survival$Study2
  K <- 2
  m <- 10
  genename <- "Simulated"
  snpnames <- sprintf("G\%s", seq(1:m))


  Fit1 <- BhGLM::bcoxph(Study1$T ~ Study1$X)
  Betah1 <- Fit1$coefficients
  Sigmah1 <- Fit1$var


  Fit2 <- BhGLM::bcoxph(Study2$T ~ Study2$X)
  Betah2 <- Fit2$coefficients
  Sigmah2 <- Fit2$var

  Betah <- list(Betah1, Betah2)
  Sigmah <- list(Sigmah1, Sigmah2)

  RES1 <- HS(Betah, Sigmah,
    kappa0 = c(0.5, 0.3), kappastar0 = c(0.5, 0.3), sigma20 = c(2, 1), s20 = c(1, 2),
    m = m, K = K, niter = 2000, burnin = 1000, nthin = 2, nchains = 2,
    a1 = 0.1, a2 = 0.1, d1 = 0.1, d2 = 0.1, c1 = 1, c2 = 1, e2 = 1, snpnames, genename
  )
  print(RES1)
}
}
\references{
Baghfalaki, T., Sugier, P. E., Truong, T., Pettitt, A. N., Mengersen, K., & Liquet, B. (2021). Bayesian meta analysis models for cross cancer genomic investigation of pleiotropic effects using group structure. Statistics in Medicine, 40(6), 1498-1518.
}
\author{
Taban Baghfalaki.
}
