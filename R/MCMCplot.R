#'  MCMC plot
#'
#'
#' @description
#'  Trace plot, density plot and ACF plot for the output of CS/DS/HS. The plot is able to draw at most ten SNPs.
#'
#'
#' @details
#'  Trace plot, density plot and ACF plot for the output of CS/DS/HS for checking convergence of MCMC chains.
#'
#'
#' @param Result All the generated results by CS/DS/HS function.
#' @param k The number of study for drawing plots, k=1,2,...,K.
#' @param nchains Number of Markov chains run in Result.
#' @param whichsnps The name of SNPs.
#'
#' @importFrom graphics lines par
#' @importFrom stats acf density
#'
#' @author Taban Baghfalaki.
#'
#' @references
#' Baghfalaki, T., Sugier, P. E., Truong, T., Pettitt, A. N., Mengersen, K., & Liquet, B. (2021). Bayesian meta analysis models for cross cancer genomic investigation of pleiotropic effects using group structure. Statistics in Medicine, 40(6), 1498-1518.
#'
#' @example inst/exampleMCMCplot.R
#'
#' @md
#' @export

MCMCplot <- function(Result = Result, k = k, nchains = nchains, whichsnps = whichsnps) {
  jjj <- k
  RES1 <- Result
  whichones <- whichsnps
  snpnames <- Result$Criteria$`Name of SNPs`
  if (length(whichones) <= 10) {
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))
    par(mfrow = c(3, length(whichones)))

    A <- RES1$MCMCChain[[1]]$Beta[[jjj]]
    colnames(A) <- snpnames
    if (nchains == 2) {
      B <- RES1$MCMCChain[[2]]$Beta[[jjj]]
      colnames(B) <- snpnames
    }
    if (nchains > 2) {
      B <- RES1$MCMCChain[[2]]$Beta[[jjj]]
      colnames(B) <- snpnames
      C <- RES1$MCMCChain[[3]]$Beta[[jjj]]
      colnames(C) <- snpnames
    }
    for (k in 1:length(whichones)) {
      P0 <- plot(A[, whichones][, k],
        type = "l",
        main = whichones[k], ylab = paste("beta", k)
      )
      # print(P0)
      if (nchains == 2) {
        L0 <- lines(B[, whichones][, k], type = "l", col = 2)
        # print(L0)
      }
      if (nchains > 2) {
        L1 <- lines(B[, whichones][, k], type = "l", col = 2)
        # print(L1)
        L2 <- lines(C[, whichones][, k], type = "l", col = 3)
        # print(L2)
      }
    }

    for (k in 1:length(whichones)) {
      P1 <- plot(density(A[, whichones][, k]), type = "l", main = whichones[k])
      # print(P1)
    }
    for (k in 1:length(whichones)) {
      ACF <- acf(A[, whichones][, k], main = whichones[k])
      # print(ACF)
    }
  } else {
    # print("Please consider at most 10 variables/SNPs.")
  }
}
