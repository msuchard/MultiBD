#include <Rcpp.h>
#include "bbd.h"

using namespace Rcpp;

// [[Rcpp::export]]
std::vector<std::complex<double>> bbd_lt_invert_Cpp(double t, const int a0, const int b0, const std::vector<double>& lambda1, const std::vector<double>& lambda2, const std::vector<double>& mu2, const std::vector<double>& gamma, const std::vector<double>& x, const std::vector<double>& y, const int A, const int B, const int maxdepth) {
  
  const double double_PI  = 3.141592653589793238463, tol = 1e-12, AA = 20.0;
  const int nblocks = 20, dim = B+1, dimsq = (B+1)*(B+1);
  int kmax = nblocks;
  std::deque<std::vector<std::complex<double>>> ig;
  std::deque<std::vector<double>> prod_mu2, prod_lambda2;
  std::vector<std::complex<double>> res((A+1-a0)*dim), phi((A+1-a0)*dimsq);
  
  for (int a=0; a<=(A-a0); ++a) {
    prod_mu2.push_back(prod_vec_Cpp(a-a0+1,A-a0,B,mu2));
    prod_lambda2.push_back(prod_vec_Cpp(a-a0+1,A-a0,B,lambda2));
  }

  for (int w=1; w<=kmax; ++w) {
    std::complex<double> s(AA/(2*t),double_PI*w/t);
    ig.push_back(bbd_lt_Cpp(s,a0,b0,lambda1,lambda2,mu2,gamma,x,y,A,B,maxdepth,phi,prod_mu2,prod_lambda2));
  }
  
  std::complex<double> s(AA/(2*t),0);
  std::vector<std::complex<double>> psum0 = bbd_lt_Cpp(s,a0,b0,lambda1,lambda2,mu2,gamma,x,y,A,B,maxdepth,phi,prod_mu2,prod_lambda2);
  
  for (int i=0;i<=(A-a0);++i)
    for(int j=0;j<=B;++j) {
      Levin levin(tol); 
      double term = 1e16, sdiff = 1e16;
      int k = 1;
      double psum = real(psum0[i*(B+1)+j])/(2*t);
      double sk,sk1;
      while ((std::abs(sdiff) > 1e-16)||(std::abs(term)>1e-3)) {
        double sgn = (k%2 == 0) ? 1.0 : -1.0;
        term = sgn*real(ig[k-1][i*(B+1)+j])/t;
        psum += term;
        double omega = k*term;
        sk = levin.next(psum,omega,1.0);
        if (k>1) sdiff = sk - sk1;
        k++;
        sk1 = sk;
        if (k > kmax) {
          for (int w=(kmax+1); w<=(kmax+nblocks); ++w) {
            std::complex<double> s(AA/(2*t),double_PI*w/t);
            ig.push_back(bbd_lt_Cpp(s,a0,b0,lambda1,lambda2,mu2,gamma,x,y,A,B,maxdepth,phi,prod_mu2,prod_lambda2));
          }
          kmax += nblocks;
        }
      }
      res[i*(B+1)+j] = sk1*exp(AA/2);
    }
  
  return(res);
}


