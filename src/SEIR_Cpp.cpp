#include "bbd.h"

// [[Rcpp::export]]
std::vector<double> SEIR_Cpp(const double t, const double alpha, const double beta, const double kappa,
                            const long int S0, const long int E0, const long int I0,
                            const int Ap1, const int Bp1, const int Cp1, const int direction,
                            const int nblocks, const double tol, int& Lmax, const int computeMode, const int nThreads) {

  const int matsize = Ap1*Bp1*Cp1;
  std::vector<double> lambda1(matsize), lambda2(matsize), lambda3(matsize);

  for (int a=0; a<Ap1; ++a) { //nSE
    for (int b=0; b<Bp1; ++b) { // nEI
      for (int c=0; c<Cp1; ++c) { // nIR
        double Spop = std::max(0.0, S0-a+0.0);
        double Epop = std::max(0.0, E0+a-b+0.0);
        double Ipop = std::max(0.0, I0+b-c+0.0);
        lambda1[a + b*Ap1 + c*Ap1*Bp1] = beta*Spop*Ipop; // Infection rate is beta*S*I
        lambda2[a + b*Ap1 + c*Ap1*Bp1] = kappa*Epop;
        lambda3[a + b*Ap1 + c*Ap1*Bp1] = alpha*Ipop;
      }
    }
  }

 return(tb_lt_invert_Cpp(t, lambda1, lambda2, lambda3, Ap1, Bp1, Cp1, direction,
                         nblocks, tol, Lmax, computeMode, nThreads));
}
