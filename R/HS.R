#' Hierarchical Spike
#'
#'
#' @description
#' Run a Gibbs sampler for a multivariate Bayesian sparse group selection model with hierarchical spike prior for detecting pleiotropic effects on the traits. This function is designed for summary statistics containing estimated regression coefficients and their estimated covariance matrices.
#'
#'
#' @details
#' Run a Gibbs sampler using a hierarchical spike.
#'
#'
#' @param Betah A list containing m-dimensional vectors of the regression coefficients for K studies.
#' @param Sigmah A list containing the positive definite covariance matrices (m*m-dimensional) which is the estimated covariance matrices of K studies.
#' @param kappa0 Initial value for kappa (its dimension is equal to nchains).
#' @param kappastar0 Initial value for kappastar (its dimension is equal to nchains).
#' @param sigma20 Initial value for sigma2 (its dimension is equal to nchains).
#' @param s20 Initial value for s2 (its dimension is equal to nchains).
#' @param m Number of variables in the group.
#' @param K  Number of traits.
#' @param niter Number of iterations for the Gibbs sampler.
#' @param burnin Number of burn-in iterations.
#' @param nthin The lag of the iterations used for the posterior analysis is defined (or thinning rate).
#' @param nchains Number of Markov chains, when nchains>1, the function calculates the Gelman-Rubin convergence statistic, as modified by Brooks and Gelman (1998).
#' @param a1,a2 Hyperparameters of kappa. Default is a1=0.1 and a2=0.1.
#' @param c1,c2 Hyperparameters of kappastar. Default is c1=0.1 and c2=0.1.
#' @param d1,d2 Hyperparameters of sigma2. Default is d1=0.1 and d2=0.1.
#' @param e2 Initial value for doing Monte Carlo EM algorithm to estimate hyperparameter of s2.
#' @param snpnames Names of variables for the group.
#' @param genename Name of group.
#'
#' @importFrom stats median pnorm quantile rbeta rbinom sd
#'
#' @return
#' - mcmcchain: The list of simulation output for all parameters.
#' - Summary: Summary statistics for regression coefficients in each study.
#' - Indicator 1: A table containing m rows of binary indicators for each study, the number of studies with nonzero signal and having pleiotropic effect by credible interval (CI). The first K columns show nonzero signals, K+1 th column includes the number of studies with nonzero signal and the last column shows an indicator for having pleiotropic effect of each SNP.
#' - Indicator 2: A table containing m rows of binary indicators for each study, the number of studies with nonzero signal and having pleiotropic effect by median. The first K columns show nonzero signals, K+1 th column includes the number of studies with nonzero signal and the last column shows an indicator for having pleiotropic effect of each SNP.
#'
#'
#'
#' @author Taban Baghfalaki.
#'
#' @references
#' Baghfalaki, T., Sugier, P. E., Truong, T., Pettitt, A. N., Mengersen, K., & Liquet, B. (2021). Bayesian meta analysis models for cross cancer genomic investigation of pleiotropic effects using group structure. Statistics in Medicine, 40(6), 1498-1518.
#'
#' @example inst/exampleHS.R
#'
#' @md
#' @export


HS <- function(Betah, Sigmah, kappa0 = kappa0, kappastar0 = kappastar0, sigma20 = sigma20, s20 = s20,
               m, K, niter = 1000, burnin = 500, nthin = 2, nchains = 2, a1 = 0.1, a2 = 0.1, d1 = 0.1, d2 = 0.1, c1 = 1, c2 = 1, e2 = 1, snpnames, genename) {
  a1 <- a1
  a2 <- a2
  d1 <- d1
  d2 <- d2
  c1 <- c1
  c2 <- c2
  e2 <- e2
  Efinal <- e2_Monte_Carlo_EM(Betah,
    Sigmah = Sigmah, kappa0 = kappa0[1],
    kappastar0 = kappastar0[1], sigma20 = sigma20[1], s20 = s20[1],
    m = m, K = K, a1 = a1, a2 = a2, d1 = d1, d2 = d2, c1 = c1, c2 = c2, e2 = e2, snpnames, genename
  )

  RES1 <- list()
  Result <- list()
  for (j in 1:nchains) {
    RES1[[j]] <- HS0(
      Betah = Betah, Sigmah = Sigmah, kappa0 = kappa0[j],
      kappastar0 = kappastar0[j], sigma20 = sigma20[j], s20 = s20[j],
      m = m, K = K, niter = niter, burnin = burnin, nthin = nthin, a1 = a1, a2 = a2, d1 = d1, d2 = d2, c1 = c1, c2 = c2,
      e2 = Efinal, snpnames = snpnames, genename = genename
    )


    Result[[j]] <- RES1[[j]]$mcmcchain
  }
  Criteria <- RES1[[j]]$Criteria
  AAT <- t(apply(RES1[[1]]$mcmcchain$Beta[[1]], 2, function(x) quantile(x, c(.025, 0.5, .975))))
  BBT <- c()
  for (k in 1:K) {
    BBT[k] <- max(RES1[[1]]$mcmcchain$Beta[[k]])
  }
  if ((nchains == 1)) {
    #  if((nchains==1)| (max(AAT)==0)| (min(BBT)==0) ){

    Summary <- list()
    for (k in 1:K) {
      Tab <- data.frame(
        snpnames, apply(RES1[[1]]$mcmcchain$Beta[[k]], 2, mean), apply(RES1[[1]]$mcmcchain$Beta[[k]], 2, sd),
        t(apply(RES1[[1]]$mcmcchain$Beta[[k]], 2, function(x) quantile(x, c(.025, 0.5, .975))))
      )
      colnames(Tab) <- cbind("Name of SNP", "Mean", "SD", "val2.5pc", "Median", "val97.5pc")
      Summary$Beta[[k]] <- Tab
    }
  } else {
    ts <- sample(1:nchains, 2)
    Summary <- list()
    for (k in 1:K) {
      AA <- list(RES1[[ts[1]]]$mcmcchain$Beta[[k]], RES1[[ts[2]]]$mcmcchain$Beta[[k]]) # Study k

      Tab <- data.frame(
        snpnames, apply(RES1[[1]]$mcmcchain$Beta[[k]], 2, mean), apply(RES1[[1]]$mcmcchain$Beta[[k]], 2, sd),
        t(apply(RES1[[1]]$mcmcchain$Beta[[k]], 2, function(x) quantile(x, c(.025, 0.5, .975)))),
        wiqid::simpleRhat(postpack::post_convert(AA))
      )
      colnames(Tab) <- cbind("Name of SNP", "Mean", "SD", "val2.5pc", "Median", "val97.5pc", "BGR")
      Summary$Beta[[k]] <- Tab
    }
  }
  RES1new <- list(MCMCChain = Result, Criteria = Criteria, Summary = Summary, Indicator = RES1[[1]]$Indicator)

  return(RES1new)
}



#' Internal: e2_Monte_Carlo_EM
#'
#' @param Betah A matrix of dimension K*m represents the regression coefficients. Each row of this matrix includes the regression coefficients for each trait.
#' @param Sigmah A symmetric block-diagonal matrix of dimension K*m is used. Each block of this matrix shows a positive definite covariance matrix which is an estimated covariance matrix of each trait.
#' @param kappa0 Initial value for kappa.
#' @param kappastar0 Initial value for kappastar.
#' @param sigma20 Initial value for sigma2.
#' @param s20 Initial value for s2.
#' @param m Number of variables in the group.
#' @param K  Number of traits.
#' @param a1,a2 Hyperparameters of kappa. Default is a1=0.1 and a2=0.1.
#' @param c1,c2 Hyperparameters of kappastar. Default is c1=0.1 and c2=0.1.
#' @param d1,d2 Hyperparameters of sigma2. Default is d1=0.1 and d2=0.1.
#' @param e2 Initial value for doing Monte Carlo EM algorithm to estimate hyperparameter of s2.
#' @param snpnames Names of variables for the group.
#' @param genename Name of group.
#'
#' @importFrom stats median pnorm quantile rbeta rbinom sd
#'

e2_Monte_Carlo_EM <- function(Betah, Sigmah, kappa0 = kappa0,
                              kappastar0 = kappastar0, sigma20 = sigma20, s20 = s20,
                              m, K, a1 = a1, a2 = a2, d1 = d1, d2 = d2, c1 = c1, c2 = c2, e2 = e2, snpnames, genename) {
  N00 <- 50
  Beta <- b <- tau <- list()
  for (k in 1:K) {
    Beta$Beta[[k]] <- b$b[[k]] <- tau$tau[[k]] <- matrix(0, N00, m)
  }

  for (k in 1:K) {
    b$b[[k]][1, ] <- Betah[[k]]
    tau$tau[[k]][1, ] <- rep(1, m)
  }

  kappa <- rep(0, N00)
  kappa[1] <- kappa0
  kappastar <- rep(0, N00)
  kappastar[1] <- kappastar0

  sigma2 <- rep(1, N00)
  sigma2[1] <- sigma20
  s2 <- rep(1, N00)
  s2[1] <- s20


  Geneplotci <- Geneplotmed <- c()

  PROB1j <- PROB2j <- c()
  PROB1 <- PROB2 <- c()

  for (r in 2:N00) {
    for (k in 1:K) {
      Vk <- diag(tau$tau[[k]][r - 1, ])
      Sigmahk <- Sigmah[[k]]

      ##################### Beta_k ####
      Omega1 <- Vk %*% MASS::ginv(Sigmahk) %*% Vk + (1 / sigma2[r - 1]) * diag(m)
      Mean1 <- MASS::ginv(Omega1) %*% Vk %*% MASS::ginv(Sigmahk) %*% matrix(Betah[[k]], m, 1)
      kappat1 <- (1 + (kappa[r - 1] / (1 - kappa[r - 1]) *
        exp(.5 * t(Mean1) %*% Omega1 %*% Mean1 -
          .5 * (determinant(Omega1, logarithm = TRUE)$modulus) - m / 2 * log(sigma2[r - 1]))))^-1

      Tab <- rbinom(1, 1, as.numeric(kappat1))
      HH1 <- MASS::ginv(Omega1)
      gdata::lowerTriangle(HH1) <- gdata::upperTriangle(HH1, byrow = TRUE)
      isSymmetric(HH1)
      b$b[[k]][r, ] <- mvtnorm::rmvnorm(1, mean = Mean1, sigma = HH1) * (1 - Tab)
      Beta$Beta[[k]][r, ] <- b$b[[k]][r, ] %*% Vk
    }

    ##################### tau1 and tau2 ####
    m1 <- v1 <- ZZ1 <- matrix(0, K, m)

    for (k in 1:K) {
      sigmabar1 <- sigmabar2 <- c()
      for (j in 1:m) {
        Sigmahk <- Sigmah[[k]]

        sigmabar1[j] <- Sigmahk[j, j] - Sigmahk[j, -j] %*% solve(Sigmahk[-j, -j]) %*% Sigmahk[-j, j]
        vkj1 <- b$b[[k]][r, j]^2 / sigmabar1[j] + 1 / s2[r - 1]
        A1 <- Betah[[k]][j] - Sigmahk[j, -j] %*% solve(Sigmahk[-j, -j]) %*% (Betah[[k]][-j] - Beta$Beta[[k]][r, -j])
        B1 <- A1 * b$b[[k]][r, j] / (vkj1 * sigmabar1[j])
        vkj1B12 <- A1^2 * b$b[[k]][r, j]^2 / (vkj1 * sigmabar1[j]^2)
        T1 <- (log(2) + pnorm(B1 * sqrt(vkj1), log.p = TRUE))
        prob1 <- (1 + kappastar[r - 1] / (1 - kappastar[r - 1]) * 1 / sqrt(vkj1 * s2[r - 1]) * exp(vkj1B12 / 2 + T1))^-1
        m1[k, j] <- B1
        v1[k, j] <- 1 / vkj1
        ZZ1[k, j] <- rbinom(1, 1, prob1)
      }

      tau$tau[[k]][r, ] <- (1 - ZZ1[k, ]) * truncnorm::rtruncnorm(m, a = 0, b = Inf, mean = m1[k, ], sd = sqrt(v1[k, ]))
    }

    ##################### kappa ####
    K0 <- 0
    for (k in 1:K) {
      if (max(b$b[[k]][r, ]) == 0) (K0 <- K0 + 1)
    }
    kappa[r] <- rbeta(1, K - K0 + a1, K0 + a2)


    ##################### kappas ####
    Tai <- c()
    for (k in 1:K) {
      Tai <- append(Tai, as.numeric(tau$tau[[k]][r, ]))
    }
    index <- 1:length(Tai)
    L0 <- length(index[Tai == 0])
    kappastar[r] <- rbeta(1, K * m - L0 + c1, L0 + c2)
    ##################### sigma2 ####
    tBB <- c()
    for (k in 1:K) {
      tBB[k] <- t(Beta$Beta[[k]][r, ]) %*% Beta$Beta[[k]][r, ]
    }
    sigma2[r] <- invgamma::rinvgamma(1, shape = m * (K - K0) / 2 + d1, rate = sum(tBB) / 2 + d2)
    while (sigma2[r] > 10) {
      sigma2[r] <- invgamma::rinvgamma(1, shape = m * (K - K0) / 2 + d1, rate = sum(tBB) / 2 + d2)
    }
    ##################### s2  ####
    tBBtau <- c()
    for (k in 1:K) {
      tBBtau[k] <- t(tau$tau[[k]][r, ]) %*% tau$tau[[k]][r, ]
    }
    s2[r] <- invgamma::rinvgamma(1, shape = (K * m - L0) / 2 + 1, rate = sum(tBBtau) / 2 + e2)
    e2 <- 1 / mean(1 / s2[1:r])
  }
  e2
}



HS0 <- function(Betah, Sigmah, kappa0, kappastar0, sigma20, s20 = s20, m, K = K, niter = 1000, burnin = 500, nthin = 2, a1 = a1, a2 = a2, c1 = c1, c2 = c2, d1 = d1, d2 = d2, e2 = e2, snpnames, genename) {
  Beta <- b <- tau <- list()
  for (k in 1:K) {
    Beta$Beta[[k]] <- b$b[[k]] <- tau$tau[[k]] <- matrix(0, niter, m)
  }

  for (k in 1:K) {
    b$b[[k]][1, ] <- Betah[[k]]
    tau$tau[[k]][1, ] <- rep(1, m)
  }

  kappa <- rep(0, niter)
  kappa[1] <- kappa0
  kappastar <- rep(0, niter)
  kappastar[1] <- kappastar0

  sigma2 <- rep(1, niter)
  sigma2[1] <- sigma20
  s2 <- rep(1, niter)
  s2[1] <- s20


  Geneplotci <- Geneplotmed <- c()

  PROB1j <- PROB2j <- c()
  PROB1 <- PROB2 <- c()

  for (r in 2:niter) {
    for (k in 1:K) {
      Vk <- diag(tau$tau[[k]][r - 1, ])
      Sigmahk <- Sigmah[[k]]

      ##################### Beta_k ####
      Omega1 <- Vk %*% MASS::ginv(Sigmahk) %*% Vk + (1 / sigma2[r - 1]) * diag(m)
      Mean1 <- MASS::ginv(Omega1) %*% Vk %*% MASS::ginv(Sigmahk) %*% matrix(Betah[[k]], m, 1)
      kappat1 <- (1 + (kappa[r - 1] / (1 - kappa[r - 1]) *
        exp(.5 * t(Mean1) %*% Omega1 %*% Mean1 -
          .5 * (determinant(Omega1, logarithm = TRUE)$modulus) - m / 2 * log(sigma2[r - 1]))))^-1
      Tab <- rbinom(1, 1, as.numeric(kappat1))
      HH1 <- MASS::ginv(Omega1)
      gdata::lowerTriangle(HH1) <- gdata::upperTriangle(HH1, byrow = TRUE)
      isSymmetric(HH1)
      b$b[[k]][r, ] <- mvtnorm::rmvnorm(1, mean = Mean1, sigma = HH1) * (1 - Tab)
      Beta$Beta[[k]][r, ] <- b$b[[k]][r, ] %*% Vk
    }

    ##################### tau1 and tau2 ####
    m1 <- v1 <- ZZ1 <- matrix(0, K, m)

    for (k in 1:K) {
      sigmabar1 <- sigmabar2 <- c()
      for (j in 1:m) {
        Sigmahk <- Sigmah[[k]]

        sigmabar1[j] <- Sigmahk[j, j] - Sigmahk[j, -j] %*% solve(Sigmahk[-j, -j]) %*% Sigmahk[-j, j]
        vkj1 <- b$b[[k]][r, j]^2 / sigmabar1[j] + 1 / s2[r - 1]
        A1 <- Betah[[k]][j] - Sigmahk[j, -j] %*% solve(Sigmahk[-j, -j]) %*% (Betah[[k]][-j] - Beta$Beta[[k]][r, -j])
        B1 <- A1 * b$b[[k]][r, j] / (vkj1 * sigmabar1[j])
        vkj1B12 <- A1^2 * b$b[[k]][r, j]^2 / (vkj1 * sigmabar1[j]^2)
        T1 <- (log(2) + pnorm(B1 * sqrt(vkj1), log.p = TRUE))
        prob1 <- (1 + kappastar[r - 1] / (1 - kappastar[r - 1]) * 1 / sqrt(vkj1 * s2[r - 1]) * exp(vkj1B12 / 2 + T1))^-1
        m1[k, j] <- B1
        v1[k, j] <- 1 / vkj1
        ZZ1[k, j] <- rbinom(1, 1, prob1)
      }

      tau$tau[[k]][r, ] <- (1 - ZZ1[k, ]) * truncnorm::rtruncnorm(m, a = 0, b = Inf, mean = m1[k, ], sd = sqrt(v1[k, ]))
    }

    ##################### kappa ####
    K0 <- 0
    for (k in 1:K) {
      if (max(b$b[[k]][r, ]) == 0) (K0 <- K0 + 1)
    }
    kappa[r] <- rbeta(1, K - K0 + a1, K0 + a2)


    ##################### kappas ####
    Tai <- c()
    for (k in 1:K) {
      Tai <- append(Tai, as.numeric(tau$tau[[k]][r, ]))
    }
    index <- 1:length(Tai)
    L0 <- length(index[Tai == 0])
    kappastar[r] <- rbeta(1, K * m - L0 + c1, L0 + c2)
    ##################### sigma2 ####
    tBB <- c()
    for (k in 1:K) {
      tBB[k] <- t(Beta$Beta[[k]][r, ]) %*% Beta$Beta[[k]][r, ]
    }
    sigma2[r] <- invgamma::rinvgamma(1, shape = m * (K - K0) / 2 + d1, rate = sum(tBB) / 2 + d2)
    while (sigma2[r] > 10) {
      sigma2[r] <- invgamma::rinvgamma(1, shape = m * (K - K0) / 2 + d1, rate = sum(tBB) / 2 + d2)
    }
    ##################### s2  ####
    tBBtau <- c()
    for (k in 1:K) {
      tBBtau[k] <- t(tau$tau[[k]][r, ]) %*% tau$tau[[k]][r, ]
    }
    s2[r] <- invgamma::rinvgamma(1, shape = (K * m - L0) / 2 + 1, rate = sum(tBBtau) / 2 + e2)
  }
  # $$$$$$$$$$$$$$$$$$$$$$$$$ New Posterior samples After Burn-in $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  N01 <- burnin + 1
  index <- N01:niter
  indexn <- index[(index %% nthin) == 0]

  kappastar <- kappastar[indexn]
  kappa <- kappa[indexn]
  sigma2 <- sigma2[indexn]
  s2 <- s2[indexn]

  for (k in 1:K) {
    Beta$Beta[[k]] <- Beta$Beta[[k]][indexn, ]
    tau$tau[[k]] <- tau$tau[[k]][indexn, ]
  }
  Beta0 <- Beta$Beta
  PGENE <- matrix(0, m, K)
  B1CI <- B2CI <- c()
  for (k in 1:K) {
    for (int in 1:m) {
      B1CI <- quantile(Beta$Beta[[k]][, int], c(0.025, 0.975))
      if ((0 < B1CI[1]) | (0 > B1CI[2])) (PGENE[int, k] <- 1)
    }
  }
  Geneplotci <- apply(PGENE, 1, sum)
  pe <- rep("No", m)
  total <- apply(PGENE, 1, sum)
  pe[total > 1] <- "Yes"
  Geneplotci <- data.frame(PGENE, total, pe)
  rownames(Geneplotci) <- snpnames
  colnames(Geneplotci) <- c(paste("Study", 1:K), "Total", "Pleiotropic effect")





  PGENE <- matrix(0, m, K)
  B1CI <- B2CI <- c()
  for (k in 1:K) {
    for (int in 1:m) {
      B1CI <- median(Beta$Beta[[k]][, int])
      if (B1CI != 0) (PGENE[int, k] <- 1)
    }
  }

  pe <- rep("No", m)
  total <- apply(PGENE, 1, sum)
  pe[total > 1] <- "Yes"
  Geneplotmed <- data.frame(PGENE, total, pe)
  rownames(Geneplotmed) <- snpnames
  colnames(Geneplotmed) <- c(paste("Study", 1:K), "Total", "Pleiotropic effect")

  mcmcchain <- list(Beta = Beta0, kappa = kappa, kappastar = kappastar, sigma2 = sigma2, s2 = s2)

  Reslast <- list("Name of Gene" = genename, "Name of SNPs" = snpnames)


  Indicator <- list("Significant studies and Variable pleiotropic effect based on CI" = Geneplotci, "Significant studies and Variable pleiotropic effect based on median" = Geneplotmed)


  OUT <- list(mcmcchain = mcmcchain, Criteria = Reslast, Indicator = Indicator)
  OUT
}
