// // [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadillo.h>
// #include "ar_coef.h"
// using namespace Rcpp;
// 
// // arma::vec ar_coef(const arma::vec, const arma::vec, const double, const int);
// // arma::vec ar_coef(arma::vec, arma::vec, double, int);
// 
// // [[Rcpp::export]]
// double sfarima_rss(const arma::vec theta,
//                    const arma::mat R_mat,
//                    const List model_order,
//                    const int K1,
//                    const int K2)
// {
//   int nar1 = as<NumericVector>(model_order["ar"])(0);
//   int nar2 = as<NumericVector>(model_order["ar"])(1);
//   int nma1 = as<NumericVector>(model_order["ma"])(0);
//   int nma2 = as<NumericVector>(model_order["ma"])(1);
//   int narma1 = nar1 + nma1;
//   int narma2 = nar2 + nma2;
// 
//   double d1 = theta(0);
//   double d2 = theta(narma1 + 1);
//   arma::colvec phi1, psi1, phi2, psi2;
//   if (nar1 > 0) {
//     phi1 = theta.subvec(1, nar1);
//   }
//   else {
//     phi1.zeros(1);
//   }
//   if (nma1 > 0) {
//     psi1 = theta.subvec(nar1 + 1, narma1);
//   }
//   else {
//     psi1.zeros(1);
//   }
//   if (nar2 > 0) {
//     phi2 = theta.subvec(narma1 + 2, narma1 + 1 + nar2);
//   }
//   else {
//     phi2.zeros(1);
//   }
//   if (nma2 > 0) {
//     psi2 = theta.subvec(narma1 + nar2 + 2, narma1 + narma2 + 1);
//   }
//   else {
//     psi2.zeros(1);
//   }
//   int n1 = R_mat.n_rows;
//   int n2 = R_mat.n_cols;
//   int Kx = std::min(K1, n1);
//   int Kt = std::min(K2, n2);
// 
//   // Function arcoef("arcoef");
//   // arma::rowvec bk1 = as<NumericVector>(arcoef(_["ar"] = phi1, _["ma"] = psi1, _["d"] = d1, _["k"] = Kx));
//   // arma::rowvec bk2 = as<NumericVector>(arcoef(_["ar"] = phi2, _["ma"] = psi2, _["d"] = d2, _["k"] = Kt));
// 
//   arma::rowvec bk1 = ar_coef(phi1, psi1, d1, Kx).t();
//   arma::rowvec bk2 = ar_coef(phi2, psi2, d2, Kt).t();
// 
//   // calculating the RSS
//   arma::mat e_est(n1, n2);
//   arma::mat xi_est(n2, n1);
//   // 1st stage
//   arma::mat tR_mat = R_mat.t();
//   for (int j = 0; j < n2; j++) {
//     if (j < Kt) {
//       xi_est.row(j) = bk2.subvec(0, j) * reverse(tR_mat.rows(0, j));
//     }
//     else {
//       xi_est.row(j) = bk2 * reverse(tR_mat.rows((j - Kt), j));
//     }
//   }
//   //2nd stage
//   arma::mat txi_est = xi_est.t();
//   for (int i = 0; i < n1; i++) {
//     if (i < Kx) {
//       e_est.row(i) = bk1.subvec(0, i) * reverse(txi_est.rows(0, i));
//     }
//     else {
//       e_est.row(i) = bk1 * reverse(txi_est.rows((i - Kx), i));
//     }
//   }
// 
//   double rss = arma::accu(arma::pow(e_est, 2));
//   return rss;
// }
// 
// 
// // /*** R
// // set.seed(1)
// // x = rnorm(10000)
// // sfarima_rss(theta = c(0.2, 0.2,  0.2,  0.3), R_mat = matrix(x, 500, 20), model_order = list(ar = c(1, 0), ma = c(1, 0)),
// //             K1 = 150, K2 = 150)
// // sfarima.rss(theta = c(0.2, 0.2, 0.2, 0.3), R_mat = matrix(x, 500, 20), model_order = list(ar = c(1, 0), ma = c(1, 0)),
// //             K1 = 150, K2 = 150)
// // */
