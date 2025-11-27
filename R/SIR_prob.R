###################################################
### Compute transition probabilities of SIR models
###################################################

#' Transition probabilities of an SIR process
#'
#' Computes the transition pobabilities of an SIR process
#' using the bivariate birth process representation
#' @param t time
#' @param alpha removal rate
#' @param beta infection rate
#' @param S0 initial susceptible population
#' @param I0 initial infectious population
#' @param nSI number of infection events
#' @param nIR number of removal events
#' @param direction direction of the transition probabilities (either \code{Forward} or \code{Backward})
#' @param power the power of the general SIR model (see Note for more details)
#' @param nblocks number of blocks
#' @param tol tolerance
#' @param computeMode computation mode
#' @param nThreads number of threads
#' @return a matrix of the transition probabilities
#' @note The infection rate and the removal rate of a general SIR model
#' are \code{beta*S^{powS}*I^{powI_inf}} and \code{alpha*I^{powI_rem}} respectively.
#' The parameter \code{power} is a list of \code{powS, powI_inf, powI_rem}.
#' Their default values are \code{powS = powI_inf = pwoI_rem = 1}, which correspond to the classic SIR model.
#' @examples
#' data(Eyam)
#'
#' loglik_sir <- function(param, data) {
#'   alpha <- exp(param[1]) # Rates must be non-negative
#'   beta  <- exp(param[2])
#'
#'   if(length(unique(rowSums(data[, c("S", "I", "R")]))) > 1) {
#'     stop ("Please make sure the data conform with a closed population")
#'   }
#'
#'   sum(sapply(1:(nrow(data) - 1), # Sum across all time steps k
#'              function(k) {
#'                log(
#'                  SIR_prob(  # Compute the forward transition probability matrix
#'                    t  = data$time[k + 1] - data$time[k], # Time increment
#'                    alpha = alpha, beta = beta,
#'                    S0 = data$S[k], I0 = data$I[k],       # From: R(t_k), I(t_k)
#'                    nSI = data$S[k] - data$S[k + 1], nIR = data$R[k + 1] - data$R[k],
#'                    computeMode = 4, nblocks = 80         # Compute using 4 threads
#'                  )[data$S[k] - data$S[k + 1] + 1,
#'                    data$R[k + 1] - data$R[k] + 1]        # To: R(t_(k+1)), I(t_(k+1))
#'                )
#'              }))
#' }
#'
#' loglik_sir(log(c(3.204, 0.019)), Eyam) # Evaluate at mode
#' @export

SIR_prob <- function(t, alpha, beta, S0, I0, nSI, nIR, direction = c("Forward","Backward"),
                     power = NULL, nblocks = 20, tol = 1e-12, computeMode = 0, nThreads = 4) {

  direction <- match.arg(direction)
  dir = 0
  if (direction == "Backward") dir = 1
  if (is.null(power)) {
    powS = 1
    powI_inf = 1
    powI_rem =1
  } else {
    powS = power$powS
    powI_inf = power$powI_inf
    powI_rem = power$powI_rem
  }

  ################################
  ### t is too small
  ### set probability 1 at a0, b0
  ################################

  if (t < tol) {
    res = matrix(0, nrow = nSI + 1, ncol = nIR + 1)
    res[1,1] = 1
    rownames(res) = 0:nSI # Infection events
    colnames(res) = 0:nIR # Removal events
    return(res)
  }

  Lmax = nblocks # initialize Lmax
  res = matrix(SIR_Cpp(t, alpha, beta, S0, I0, nSI + 1, nIR + 1, dir,
                       powS, powI_inf, powI_rem,
                       nblocks, tol, Lmax, computeMode, nThreads),
               nrow = nSI + 1, byrow = T)

  rownames(res) = 0:nSI # Infection events
  colnames(res) = 0:nIR # Removal events

  return(abs(res))
}
