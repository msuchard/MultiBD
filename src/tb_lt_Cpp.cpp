#include "bbd.h"

void tb_lt_Cpp(const mytype::ComplexNumber s, const std::vector<double>& lambda1,
               const std::vector<double>& lambda2, const std::vector<double>& lambda3,
               const int Ap1, const int Bp1, const int Cp1, const int direction,
               const std::vector<double>& yvec, mytype::ComplexVector& f) {

  /////////////////////////////////////////////////////////////////////////////////
  // The following code computes Laplace transform of the transition probabilities
  /////////////////////////////////////////////////////////////////////////////////

  if (direction == 0) { // Forward direction
    // for (int i=0; i<Ap1; ++i) {
    //   for (int j=0; j<Bp1; ++j) {
    //     for (int k=0; k<Cp1; ++k) {
    //       f[i + j*Ap1 + k*Ap1*Bp1] = zero;
    //     }
    //   }
    // }
    f[0] = 1; // Assume that all other elements of f is 0
    for (int i=0; i<Ap1; ++i) {
      for (int j=0; j<Bp1; ++j) {
        for (int k=0; k<Cp1; ++k) {
          int tmp = i + j*Ap1 + k*Ap1*Bp1;
          f[tmp] = f[tmp]/(s + yvec[tmp]);
          if (i < (Ap1-1)) f[tmp + 1] = f[tmp + 1] + lambda1[tmp]*f[tmp];
          if (j < (Bp1-1)) f[tmp + Ap1] = f[tmp + Ap1] + lambda2[tmp]*f[tmp];
          if (k < (Cp1-1)) f[tmp + Ap1*Bp1] = f[tmp + Ap1*Bp1] + lambda3[tmp]*f[tmp];
        }
      }
    }
  } else { // Backward direction
    f[Ap1*Bp1*Cp1 - 1] = 1; // Assume that all other elements of f is 0
    for (int i=(Ap1-1); i>=0; --i) {
      for (int j=(Bp1-1); j>=0; --j) {
        for (int k=(Cp1-1); k>=0; --k) {
          int tmp = i + j*Ap1 + k*Ap1*Bp1;
          f[tmp] = f[tmp]/(s + yvec[tmp]);
          if (i > 0) f[tmp - 1] = f[tmp - 1] + lambda1[tmp - 1]*f[tmp];
          if (j > 0) f[tmp - Ap1] = f[tmp - Ap1] + lambda2[tmp - Ap1]*f[tmp];
          if (k > 0) f[tmp - Ap1*Bp1] = f[tmp - Ap1*Bp1] + lambda3[tmp - Ap1*Bp1]*f[tmp];
        }
      }
    }
  }

}
