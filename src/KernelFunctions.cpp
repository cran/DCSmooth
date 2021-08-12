// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "DCSmooth_types.h"

using namespace Rcpp;

//-----------------------Weight Functions for LPR-----------------------------//

// Truncated
arma::vec weights_T(arma::vec& u, double q = 1, int mu = 2)
{
  arma::vec weights_out{ pow(1 + u, mu) % pow(1 - u, mu) };
  return weights_out;
}

// Müller
arma::vec weights_M(arma::vec& u, double q = 1, int mu = 2)
{
  arma::vec weights_out{ pow(1 + u, mu) % pow(q - u, mu) };
  return weights_out;
}

// Müller-Wang
arma::vec weights_MW(arma::vec& u, double q = 1, int mu = 2)
{
  arma::vec weights_out{ pow(1 + u, mu) % pow(q - u, std::max(mu - 1, 0)) };
  return weights_out;
}

// [[Rcpp::export]]
XPtr<weightPtr> weight_fcn_assign(std::string fstr) {
  if (fstr == "T")
    return(XPtr<weightPtr>(new weightPtr(&weights_T)));
  else if (fstr == "M")
    return(XPtr<weightPtr>(new weightPtr(&weights_M)));
  else if (fstr == "MW")
    return(XPtr<weightPtr>(new weightPtr(&weights_MW)));
  else
    return XPtr<weightPtr>(R_NilValue); // runtime error as NULL no XPtr
}

//----------------------------Kernel Functions--------------------------------//

// Müller-Wang kernels
// Form is: Type_k_mu_nu

arma::vec kernFkt_MW200(arma::vec& u, double q = 1)
{
  arma::vec qVec{ arma::vec(u.n_rows).fill(q) };
  arma::vec q2{ pow(qVec, 2) };
  double n0{ 2.0/pow(q + 1, 3) };

  return n0 * (2*(q2 - qVec + 1) + 3*u%(1 - qVec));
}

// [[Rcpp::export]]
arma::vec kernFkt_MW210(arma::vec& u, double q = 1)
{
  arma::vec qVec{ arma::vec(u.n_rows).fill(q) };
  arma::vec wq{ (1 + u) };
  arma::vec q2{ pow(qVec, 2) };
  double n0{ 6.0/pow(q + 1, 4) };

  return n0 * ((3*q2 - 2*qVec + 1) + 2*u % (1 - 2*qVec)) % wq;
}

// [[Rcpp::export]]
arma::vec kernFkt_MW220(arma::vec& uVec, double q)
{
  arma::uword nBound{ uVec.size() };
  
  arma::vec uOut(nBound);
  
  double q2{ q * q };
  double q2_1{ pow(1 + q, 6) };
  double out{ };
  double u{ };
  
  for (int i{ 0 }; i < nBound; ++i)
  {
    u   = uVec(i);
    out = (60.0/q2_1) * (u*(2 - 3*q) + (1 - 2*q + 2*q2))
                      * (1 + u)*(1 + u)*(q - u);
    uOut(i) = out;
  }
  return uOut;
}

// [[Rcpp::export]]
arma::vec kernFkt_MW320(arma::vec& u, double q = 1)
{
  arma::vec qVec{ arma::vec(u.n_rows).fill(q) };
  arma::vec wq{ (1 + u) % (1 + u) % (qVec - u) };
  arma::vec q4{ pow(qVec, 4) };
  arma::vec q3{ pow(qVec, 3) };
  arma::vec q2{ pow(qVec, 2) };
  double n0{ 60.0/pow(q + 1, 8) };
  arma::vec mid{ (10*q4 - 30*q3 + 39*q2 - 16*qVec + 3)
                  - 7*(5*qVec - 2) % (qVec - 1) % (qVec - 1) % u
                  + 14*(2*q2 - 4*qVec + 1) % u % u };

  return (n0 * mid % wq);
}

// [[Rcpp::export]]
arma::vec kernFkt_MW420(arma::vec& uVec, double q)
{
  arma::uword nBound{ uVec.size() };
  
  arma::vec uOut(nBound);
  
  double q2{ q * q };
  double q3{ q2 * q };
  double q4{ q3 * q };
  double q5{ q4 * q };
  double q6{ q5 * q };
  double q10_1{ pow(1 + q, 10) };
  double out{ };
  double u{ };
  
  for (int i{ 0 }; i < nBound; ++i)
  {
    u   = uVec(i);
    out = (420.0/q10_1) * ( u*u*u*(12 - 90*q + 120*q2 - 30*q3) 
                  + u*u*(18 - 144*q + 300*q2 - 240*q3 + 54*q4)
                  + u*(8 - 70*q + 216*q2 - 280*q3 + 152*q4 - 30*q5)
                  + (1 - 10*q + 45*q2 - 84*q3 + 77*q4 - 30*q5 + 5*q6))
                * (1 + u)*(1 + u)*(q - u);
    uOut(i) = out;
  }
  return uOut;
}

// [[Rcpp::export]]
arma::vec kernFkt_MW421(arma::vec& uVec, double q)
{
  arma::uword nBound{ uVec.size() };
  
  arma::vec uOut(nBound);
  
  double q2{ q * q };
  double q3{ q2 * q };
  double q4{ q3 * q };
  double q5{ q4 * q };
  double q10_1{ pow(1 + q, 10) };
  double out{ };
  double u{ };
  double u2{ };
  double u3{ };
  
  for (int i{ 0 }; i < nBound; ++i)
  {
    u   = uVec(i);
    u2  = u * u;
    u3  = u2 * u;
    
    out = 840*(q - u)*(u + 1)*(u + 1) * ( -15*q5 +
             q4*(97*u + 76) -
             q3*(183*u2 + 344*u + 140) +
             q2*(105*u3 + 480*u2 + 444*u + 108) -
             q *(210*u3 + 381*u2 + 212*u + 35) +
                (63*u3 + 90*u2 + 37*u + 4)) / q10_1;
      uOut(i) = out;
  }
  return uOut;
}

// [[Rcpp::export]]
arma::vec kernFkt_MW422(arma::vec& uVec, double q)
{
  arma::uword nBound{ uVec.size() };
  
  arma::vec uOut(nBound);
  
  double q2{ q * q };
  double q3{ q2 * q };
  double q4{ q3 * q };
  double q10_1{ pow(1 + q, 10) };
  double out{ };
  double u{ };
  
  for (int i{ 0 }; i < nBound; ++i)
  {
    u   = uVec(i);
    out = (5040.0/q10_1) * ( u*u*u*(56 - 70*q) + u*u*(77 - 182*q + 119*q2)
                         + u*(30 - 127*q + 160*q2 - 61*q3)
                         + (3 - 24*q + 50*q2 - 40*q3 + 9*q4))
                       * (1 + u)*(1 + u)*(q - u);
    uOut(i) = out;
  }
  return uOut;
}

// Truncated Kernels
// Form is 

// [[Rcpp::export]]
arma::vec kernFkt_TR420(arma::vec& uVec, double q)
{
  arma::uword nBound{ uVec.size() };
  
  arma::vec uOut(nBound);
  
  double q2{ q * q };
  double q3{ q2 * q };
  double q4{ q3 * q };
  double q5{ q4 * q };
  double q6{ q5 * q };
  double q7{ q6 * q };
  double q8{ q7 * q };
  double q9{ q8 * q };
  double q10{ q9 * q };
  double q11{ q10 * q };
  double q12{ q11 * q };
  double q_9{ pow(q + 1, 9) };
  double out{ };
  double u{ };
  double u2{ };
  double u3{ };
  
  for (int i{ 0 }; i < nBound; ++i)
  {
    u  = uVec(i);
    u2 = u * u;
    u3 = u2* u;
    out = 840*pow(u - 1, 2)*pow(u + 1, 2) * (1960*q12 - 
              q11*(8820*u      + 29400)       +
              q10*(12600*u2   + 132300*u     + 203952) -
              q9 *(5775*u3    + 189000*u2   + 905310*u + 864080) +
              q8 *(86625*u3   + 1280160*u2  + 3701250*u + 2477160) - 
              q7 *(582120*u3  + 5090400*u2  + 9930060*u + 4982040) +
              q6 *(2263800*u3 + 12910200*u2 + 18132660*u + 7056480) -
              q5 *(5470542*u3 + 21486600*u2 + 22702176*u + 6950304) +
              q4 *(8322930*u3 + 23479920*u2 + 19202400*u + 4644000) -
              q3 *(7780080*u3 + 16440720*u2 + 10553760*u + 2012000) +
              q2 *(4158000*u3 + 6923520*u2  + 3497760*u + 514560) -
              q  *(1063755*u3 + 1500480*u2  + 606690*u + 64320) +
                  (70917*u3   + 100032*u2   + 40446*u + 4288)) / 
                  (q_9*(105*q8 - 2520*q7 + 26236*q6 - 146664*q5 +
                  479670*q4 - 941800*q3 + 1089660*q2 - 682968*q + 178537));
    uOut(i) = out;
  }
  return uOut;
}

// [[Rcpp::export]]
arma::vec kernFkt_TR422(arma::vec& uVec, double q)
{
  arma::uword nBound{ uVec.size() };
  
  arma::vec uOut(nBound);
  
  double q2{ q * q };
  double q3{ q2 * q };
  double q4{ q3 * q };
  double q5{ q4 * q };
  double q6{ q5 * q };
  double q7{ q6 * q };
  double q8{ q7 * q };
  double q9{ q8 * q };
  double q10{ q9 * q };
  double q_9{ pow(q + 1, 9) };
  double out{ };
  double u{ };
  double u2{ };
  double u3{ };
  
  for (int i{ 0 }; i < nBound; ++i)
  {
    u  = uVec(i);
    u2 = u * u;
    u3 = u2* u;
    out = 5040*pow(u - 1, 2)*pow(u + 1, 2)*(4200*q10 -
              q9*(19845*u + 63000*q9) +
              q8*(29400*u2 + 297675*u + 426720) -
              q7*(13860*u3 + 441000*u2 + 1968120*u - 1696800) +
              q6*(207900*u3 + 2865520*u2 + 7295400*u + 4303400) -
              q5*(1333332*u3 + 10054800*u2 + 16618266*u + 7162200) +
              q4*(4476780*u3 + 20784120*u2 + 24117030*u + 7826640) -
              q3*(8468460*u3 + 26072200*u2 + 22251600*u + 5480240) +
              q2*(9050580*u3 + 19568640*u2 + 12514320*u + 2307840) -
              q *(5086620*u3 + 8137920*u2 + 3848985*u + 500160) +
                 (1167012*u3 + 1460032*u2 + 482391*u + 33344)) /
                (q_9*(105*q8 - 2520*q7 + 26236*q6 - 146664*q5 +
                  479670*q4 - 941800*q3 + 1089660*q2 - 682968*q + 178537));
    uOut(i) = out;
  }
  return uOut;
}

// [[Rcpp::export]]
XPtr<funcPtr> kernel_fcn_assign(std::string fstr) {
  if (fstr == "MW_200")
    return(XPtr<funcPtr>(new funcPtr(&kernFkt_MW200)));
  else if (fstr == "MW_210")
    return(XPtr<funcPtr>(new funcPtr(&kernFkt_MW210)));
  else if (fstr == "MW_220")
    return(XPtr<funcPtr>(new funcPtr(&kernFkt_MW220)));
  else if (fstr == "MW_320")
    return(XPtr<funcPtr>(new funcPtr(&kernFkt_MW320)));
  else if (fstr == "MW_420")
    return(XPtr<funcPtr>(new funcPtr(&kernFkt_MW420)));
  else if (fstr == "MW_422")
    return(XPtr<funcPtr>(new funcPtr(&kernFkt_MW422)));
  else
    return XPtr<funcPtr>(R_NilValue); // runtime error as NULL no XPtr
}

// [[Rcpp::export]]
arma::vec kernel_fcn_use(arma::vec x, double q, SEXP xpsexp) {
  XPtr<funcPtr> xpfun(xpsexp);
  funcPtr fun = *xpfun;
  arma::vec y = fun(x, q);
  return y;
}
