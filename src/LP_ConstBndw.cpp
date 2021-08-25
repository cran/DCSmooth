// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include "DCSmooth.h"
#include "DCSmooth_types.h"

using namespace Rcpp;

//---------------------------------------------------------------------------//

// [[Rcpp::export]]
arma::mat LPSmooth_matrix2(const arma::mat yMat, const double h,
                           const int polyOrder, const int drv, SEXP kernFcnPtr)
{
  // get additional information on nX, nT, bndw etc.
  int nRow{ static_cast<int>(yMat.n_rows) };
  int nCol{ static_cast<int>(yMat.n_cols) };
  int bndw{ std::max(static_cast<int>(h * nCol), polyOrder + 1) };
                        // calculate absolute bandwidth, decimals will be dumped
  int windowWidth{ 2*bndw + 1 };  // width of estimation window

  // result matrix
  arma::mat yMatOut(nRow, nCol);

  // enable Kernel function
  XPtr<funcPtr> xpfun(kernFcnPtr);
  funcPtr kernFcn = *xpfun;

  // calculate weights
  arma::colvec  uVec{ arma::regspace(bndw, -bndw)/std::max(h * nCol + 1,
                        static_cast<double>(polyOrder + 1)) }; // vector from -1 to 1 to compute weights
  arma::colvec  xVec{ arma::regspace(-bndw, bndw)/std::max(h * nCol,
                        static_cast<double>(polyOrder + 1)) }; // vector from -1 to 1 to compute weights
  arma::colvec  weightsVec{ (kernFcn(uVec, 1)) };         // computation of weights // sqrt removed

  // smoothing over boundaries
  for (int colIndex{ 0 }; colIndex < bndw; ++colIndex)
  {
    if (colIndex == (nCol/2) + 1)     // break condition, if bndw > 0.5*nCol
    {
      break;
    }

    int startIndex{ bndw - colIndex };
    int stopIndex{ std::min(windowWidth - 1, startIndex + nCol - 1) };

    // calculate weights for linear regression
    arma::colvec weightsBound{ weightsVec.subvec(startIndex, stopIndex) };
    arma::colvec xBound{ xVec.subvec(startIndex, stopIndex) };
    arma::mat    xMatBound{ xMatrix(xBound, polyOrder) };
    arma::mat    xWeightBound{ weightMatrix(weightsBound, xMatBound) };
    arma::mat    xMatSolved{ arma::inv(xWeightBound.t() * xMatBound)
                                * xWeightBound.t() };
    arma::rowvec xWeightsLeft{ xMatSolved.row(drv) * factorialFunction(drv) /
                               std::pow(h, drv) };
    // pow(-1, drv) ensures the correct sign
    arma::rowvec xWeightsRight{ std::pow(-1, drv) * xWeightsLeft };

    // calculation of estimates (complete column)
    yMatOut.col(colIndex) = yMat.cols(0, xBound.n_rows - 1) * xWeightsLeft.t();
    yMatOut.col(nCol - colIndex - 1) = arma::reverse(yMat.cols(nCol -
                              xBound.n_rows, nCol - 1), 1) * xWeightsRight.t();
  }

  if (h < 0.5)
  {
    // smothing over interior values
    arma::mat yInterior(nRow, windowWidth); // empty vector for use inside loop
    arma::mat coefMat(polyOrder, windowWidth);    // empty matrix for lm results
    arma::mat xMatInterior{ xMatrix(xVec, polyOrder) };  // compute x-matrix for lm-regression
    arma::mat xMatWeight{ weightMatrix(weightsVec, xMatInterior) };
    arma::mat xMatSolved{ arma::inv(xMatWeight.t() * xMatInterior)
      * xMatWeight.t() };         // compute inv(X*W*X)^(-1)*X*W only once here
    arma::rowvec xWeightsVec{ xMatSolved.row(drv) *  factorialFunction(drv) /
                              std::pow(h, drv) };

    // Loops smooth over the columns, conditional on rows. That is, every row is
    // consiedered to be an individual time series. To speed up computation, the
    // smoothing order is inverted (computing of weights only once per column, as
    // the weights are the same on a grid)
    for (int colIndex{ bndw }; colIndex < (nCol - bndw); ++colIndex)                  // outer loop over columns
    {
      yInterior  = yMat.cols(colIndex - bndw, colIndex + bndw);
      yMatOut.col(colIndex) = yInterior * xWeightsVec.t();
    }
  }

  return yMatOut;
}

//---------------------------------------------------------------------------//

// [[Rcpp::export]]
arma::mat LP_dcs_const0(arma::mat yMat, arma::colvec hVec,
                       arma::icolvec polyOrderVec, arma::icolvec drvVec,
                       SEXP kernFcnPtr_x, SEXP kernFcnPtr_t)
{
  arma::mat mMatTemp{ LPSmooth_matrix2(yMat, hVec(1), polyOrderVec(1),
                                       drvVec(1), kernFcnPtr_t) };
  arma::mat yMatOut{ LPSmooth_matrix2(mMatTemp.t(), hVec(0), polyOrderVec(0),
                                      drvVec(0), kernFcnPtr_x) };

  return yMatOut.t();
}
