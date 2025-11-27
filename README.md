MultiBD
======

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/MultiBD)](http://cran.r-project.org/package=MultiBD)


`MultiBD` is an `R` package for direct likelihood-based inference of multivariate birth-death processes. 

## Installation

1. Install `CRAN` release version:
```{r}
install.packages("MultiBD")
```

2. Install the bleeding-edge version of `MultiBD` from `github`:
```{r}
devtools::install_github("msuchard/MultiBD")
```

## Short example

```{r}
library(MultiBD)
data(Eyam)

loglik_sir <- function(param, data) {
  alpha <- exp(param[1]) # Rates must be non-negative
  beta  <- exp(param[2])
  
  # Set-up SIR model
  drates1 <- function(a, b) { 0 }
  brates2 <- function(a, b) { 0 }
  drates2 <- function(a, b) { alpha * b     }
  trans12 <- function(a, b) { beta  * a * b }
  
  sum(sapply(1:(nrow(data) - 1), # Sum across all time steps k
             function(k) {
               log(
                 dbd_prob(  # Compute the transition probability matrix
                   t  = data$time[k + 1] - data$time[k], # Time increment
                   a0 = data$S[k], b0 = data$I[k],       # From: S(t_k), I(t_k)                                      
                   drates1, brates2, drates2, trans12,
                   a = data$S[k + 1], B = data$S[k] + data$I[k] - data$S[k + 1],
                   computeMode = 4, nblocks = 80         # Compute using 4 threads
                 )[1, data$I[k + 1] + 1]                 # To: S(t_(k+1)), I(t_(k+1))
               )
             }))
}

loglik_sir(log(c(3.204, 0.019)), Eyam) # Evaluate at mode
```


## Vignettes

1. [Simple MCMC under SIR](https://github.com/msuchard/MultiBD/blob/master/inst/doc/SIR-MCMC.pdf)
2. [SIR model and proposed branching approximation](https://github.com/msuchard/MultiBD/blob/master/inst/doc/SIRtrans.pdf)

## License
`MultiBD` is licensed under Apache License 2.0

## Development status

[![Build Status](https://travis-ci.org/msuchard/MultiBD.svg?branch=master)](https://travis-ci.org/msuchard/MultiBD)

Beta

## Acknowledgements
- This project is supported in part through the National Science Foundation grant DMS 1264153 and National Institutes of Health grant R01 AI107034.

## References

1. Ho LST, Xu J, Crawford FW, Minin VN, Suchard MA (2018).
[Birth/birth-death processes and their computable transition probabilities with biological applications](https://link.springer.com/article/10.1007/s00285-017-1160-3).
Journal of Mathematical Biology 76(4) 911-944.
2. Ho LST, Crawford FW, Suchard MA (2018).
[Direct likelihood-based inference for discretely observed stochastic compartmental models of infectious disease](https://arxiv.org/abs/1608.06769).
Annals of Applied Statistics. In press.
