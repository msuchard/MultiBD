## Package patch following email from Brian Ripley noting compiler deprecation of C++14

1. Modernized author-list in `DESCRIPTION`
2. Fixed C++ deprecation of `typename std::result_of<F(Args...)>::type`

## Test environments
* local macOS install (arm64), R 4.4
* rhub (multiple architectures), R 4.5
* win-builder (devel, release, oldrel)

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs.

## Downstream dependencies
There are currently no downstream dependencies.
