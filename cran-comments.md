## Package patch following email from Brian Ripley noting:
  `warning: implicit capture of ‘this’ via ‘[=]’ is deprecated in C++20`

1. Fixed C++20 deprecation of `implicit capture of this`

## Test environments
* local macOS install (arm64), R 4.4
* rhub (multiple architectures), R 4.5
* win-builder (devel, release, oldrel)

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs.

## Downstream dependencies
There are currently no downstream dependencies.
