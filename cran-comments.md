## Package patch following email from Brian Ripley noting:
  ` warning: ignoring return value of function declared with 'nodiscard' attribute`

1. Fixed C++ warning of `ignoring return value`

## Test environments
* local macOS install (arm64), R 4.4
* rhub (multiple architectures), R 4.5
* win-builder (devel, release, oldrel)

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs.

## Downstream dependencies
There are currently no downstream dependencies.
