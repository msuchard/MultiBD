###################################################
### Compute transition probabilities of SEIR models
###################################################

#' Transition probabilities of an SEIR process
#'
#' Computes the transition pobabilities of an SIR process
#' using the bivariate birth process representation
#' @param t time
#' @param alpha removal rate
#' @param kappa rate at which an exposed person becomes infective
#' @param beta infection rate
#' @param S0 initial susceptible population
#' @param E0 initial exposed population
#' @param I0 initial infectious population
#' @param nSE number of infection events
#' @param nEI number of events at which an exposed person becomes infective
#' @param nIR number of removal events
#' @param direction direction of the transition probabilities (either \code{Forward} or \code{Backward})
#' @param nblocks number of blocks
#' @param tol tolerance
#' @param computeMode computation mode
#' @param nThreads number of threads
#' @return a matrix of the transition probabilities
#' @export

SEIR_prob <- function(t, alpha, beta, kappa, S0, E0, I0, nSE, nEI, nIR, direction = c("Forward","Backward"),
                     nblocks = 20, tol = 1e-12, computeMode = 0, nThreads = 4) {

  direction <- match.arg(direction)
  dir = 0
  if (direction == "Backward") dir = 1

  ####################################
  ### t is too small
  ### set probability 1 at a0, b0, c0
  ####################################

  res = array(0, dim= c(nSE + 1,nEI + 1,nIR+1))
  if (t < tol) {
    res[1,1,1] = 1
  } else {
    Lmax = nblocks # initialize Lmax
    tmp = SEIR_Cpp(t, alpha, beta, kappa, S0, E0, I0, nSE + 1, nEI + 1, nIR + 1, dir, nblocks, tol,
                         Lmax, computeMode, nThreads)
    for (i in 0:nSE) {
      for (j in 0:nEI) {
        for (k in 0:nIR)
          res[i+1,j+1,k+1] = tmp[i + j*(nSE+1) + k*(nSE+1)*(nEI+1) + 1]
      }
    }
  }

  return(abs(res))
}


