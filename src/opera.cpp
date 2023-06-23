#include<RcppEigen.h>
#include<Rcpp.h>
#include<vector>
#include<cmath>
#include<iostream>
#include<limits>
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

//' C++ function to call the R function [quadprog::solve.QP()]
//'
//' This function serves as a bridge between C++ and R, allowing you to call the R function [quadprog::solve.QP()] from C++ code.
//'
//' @param Dmat A matrix appearing in the quadratic function to be minimized.
//' @param dvec A vector appearing in the quadratic function to be minimized.
//' @param Amat A matrix defining the constraints under which we want to minimize the quadratic function.
//'
//' @seealso [quadprog::solve.QP()]
//'
//' @return A list containing the solution obtained from the R function [quadprog::solve.QP()].
//'@export
//[[Rcpp::export]]
List quadprog_solveR(const MatrixXd & Dmat, const VectorXd & dvec,
                     const MatrixXd & Amat, const VectorXd & bvec){
  List result;
  Environment myEnv = Environment::namespace_env("quadprog");
  Function quadR = myEnv["solve.QP"];
  result = quadR(Dmat, dvec, Amat, bvec);
  return result;
}

//' Call R Function to Compute Eigenvalues and Eigenvectors
//'
//' This C++ function calls an R function to compute the eigenvalues and eigenvectors of a numeric (double, integer, logical) or complex matrix.
//'
//' @param mat An approximately positive definite matrix.
//' @return A list containing the eigenvalues and eigenvectors.
//'
//'@export
//[[Rcpp::export]]
List ED(const MatrixXd & mat){
  List result;
  Environment myEnv = Environment::namespace_env("base");
  Function EG = myEnv["eigen"];
  result = EG(mat);
  return result;
}

//' Replace Negative Values with Zeros in a Vector
//'
//' This function replaces all the negative values of a vector with zeros.
//'
//' @param v A numeric vector.
//' @return The non-negative version of the vector.
//'
//'@export
//[[Rcpp::export]]
VectorXd nNv(const VectorXd & v){
  VectorXd nv(v.size());
  for(int i = 0 ; i < v.size(); i++){
    if(v[i] >= 0){
      nv[i] = v[i];
    }else{
      nv[i] = 0;
    }
  }
  return(nv);
}


//' Compute Nearest Positive Definite Matrix
//'
//' This function computes the nearest positive definite matrix to an approximately positive definite matrix.
//'
//' @param mat An approximately positive definite matrix.
//'
//' @return The nearest positive definite matrix.
//'
//'@export
//[[Rcpp::export]]
MatrixXd CnearPD(const MatrixXd & mat){
  List r;
  MatrixXd pdm;
  r = ED(mat);
  MatrixXd ev = r.at(1);
  VectorXd nv = nNv(r.at(0));

  for(int k = 0; k < nv.size(); k++){
    if(nv[k] == 0)
      nv[k] = 0.0001;
  }

  pdm = ev*nv.asDiagonal()*ev.transpose();
  return(pdm);
}


//' Calculate Sum of Vector Elements
//'
//' This function calculates the sum of elements in a given vector.
//'
//' @param v A vector.
//' @return The sum of the vector elements.
//'
//'@export
//[[Rcpp::export]]
double Sum(const VectorXd & v){
  double s = 0;
  for(int i = 0; i < v.size(); i++){
    s+=v[i];
  }
  return(s);
}

//' Calculate Minimum of Vector Elements
//'
//' This function calculates the minimum value among the elements of a given vector.
//'
//' @param v A vector.
//' @return The minimum value of the vector elements.
//'
//'@export
//[[Rcpp::export]]
double Min(const VectorXd & v){
  double m = v[0];
  for(int i = 0; i < v.size(); i++){
    if(v[i] < m)
      m = v[i];
  }
  return(m);
}

//' Check for Inf or NaN in Vector
//'
//' This function checks whether a given vector contains any elements that are either Inf (infinity) or NaN (not-a-number).
//'
//' @param v A vector.
//' @return A boolean value indicating whether the vector contains Inf or NaN.
//'
//'@export
//[[Rcpp::export]]
bool wInf(const VectorXd & v){
  for(int i = 0; i < v.size(); i++){
    if(isinf(v[i])||(v[i] != v[i]))
      return(isinf(v[i])||(v[i] != v[i]));
  }
  return(false);
}

//' Round Vector Elements
//'
//' This function rounds the elements of a given vector to the nearest integer, based on a specified threshold.
//'
//' @param v A vector.
//' @param eps A small number representing the rounding threshold.
//' @return A vector with rounded elements.
//'
//'@export
//[[Rcpp::export]]
VectorXd Round(const VectorXd & v, double eps){
  VectorXd r(v.size());
  for(int i = 0; i < v.size(); i++){
    if(abs(v[i]) <= eps)
      r[i] = 0;
    else
      r[i] = round(v[i]/eps)*eps;
  }
  return(r);
}

//' Count Non-Zero Distinct Values in Vector
//'
//' This function calculates the number of non-zero distinct values in a given vector.
//'
//' @param v A vector.
//' @return The number of non-zero distinct values in the vector.
//'
//'@export
//[[Rcpp::export]]
int countDistinct(const VectorXd & v){
  int size = v.size();
  int count = 1;
  if(size <= 1){
    count = size;
  }else{
    for (int i = 1; i < size; i++){
      if ((v.head(i).array() == v[i]).count() == 0){
        count = count + 1;
      }
    }
  }
  if((v.array() == 0).any()){
    count = count - 1;
  }
  return(count);
}

//' Calculate Cumulative Sums of Vector Elements
//'
//' This function calculates the cumulative sums of elements in a given vector.
//'
//' @param v A vector.
//' @return A vector containing the cumulative sums of the input vector.
//'
//'@export
//[[Rcpp::export]]
VectorXd cumSum(const VectorXd  & v){

  VectorXd s(v.size());
  double is;

  for(int j = 0; j < v.size(); j++){
    is = 0;
    for(int k = 0; k < j+1; k++){
      is += v[k];
    }
    s[j] = is;
  }
  return(s);
}

//' Calculate Reverse Cumulative Sums of Vector Elements
//'
//' This function calculates the reverse cumulative sums of elements in a given vector.
//'
//' @param v A vector.
//' @return A vector containing the reverse cumulative sums of the input vector.
//'
//'@export
//[[Rcpp::export]]
VectorXd revCumSum(const VectorXd  &  v){

  VectorXd rv(v.size());
  VectorXd s(v.size());
  VectorXd rs(v.size());
  double is;

  for(int i = 0; i < v.size(); i++){
    rv[i] = v[v.size()-i-1];
  }


  for(int j = 0; j < rv.size(); j++){
    is = 0;
    for(int k = 0; k < j+1; k++){
      is += rv[k];
    }
    s[j] = is;
  }

  for(int l = 0; l < s.size(); l++){
    rs[l] = s[s.size()-l-1];
  }

  return(rs);
}

//' Create Equally-Spaced Logarithmic Sequence
//'
//' This function creates an equally-spaced sequence in the log scale, also known as a geometric series, between a starting value and an ending value.
//'
//' @param a A numeric value representing the starting value.
//' @param b A numeric value representing the ending value.
//' @param c The length of the sequence.
//' @return An equally-spaced sequence in the log scale (a geometric series).
//'
//'@export
//[[Rcpp::export]]
VectorXd gs(const double & a, const double & b, const int & c = 30){
  double d = (log(b) - log(a))/(c-1);
  VectorXd v(c);
  for(int i = 0; i < c; i++){
    v[i] = exp(d*i + log(a));
  }
  return(v);
}


//' Calculate Negative Log (Partial) Likelihood
//'
//' This function calculates the negative log (partial) likelihood based on the given parameters.
//'
//' @param Z A matrix of nodes.
//' @param cen A vector of censoring statuses or binary outcomes.
//' @param beta A vector of node coefficients.
//' @param y A vector of outcomes.
//' @param withCov Whether any covariates need adjusting.
//' @param cov A matrix of covariates.
//' @param theta A vector of covariate coefficients.
//' @param type A string indicating the type of outcome, currently either "survival" or "binary".
//'
//' @return The negative log likelihood.
//'
//'@export
//[[Rcpp::export]]
double logLK(const MatrixXd & Z, const VectorXd & cen,
             const VectorXd & beta, const VectorXd & y,
             const MatrixXd & cov, const VectorXd & theta,
             const std::string type = "surv", const bool withCov = false){

  VectorXd eta;
  VectorXd logeta;
  VectorXd expeta;
  double r;

  if(withCov)
    eta = Z*beta + cov*theta;
  else
    eta = Z*beta;

  expeta = exp(eta.array());

  if(type == "bin"){
    logeta = log( (VectorXd::Ones(expeta.size()) + expeta).array() );
    r  =  -y.transpose() * eta + Sum(logeta);
  }else{
    logeta = log(revCumSum(expeta).array());
    r = logeta.dot(cen) - cen.dot(eta);
  }
  return(r);
}

//' Calculate First and Second Derivatives of Log (Partial) Likelihood
//'
//' This function calculates the first and second derivatives of the log (partial) likelihood based on the given parameters.
//'
//' @param Z A matrix of nodes.
//' @param cen A vector of censoring statuses or binary outcomes.
//' @param beta A vector of node coefficients.
//' @param y A vector of outcomes.
//' @param cov Covariates.
//' @param theta Covariate coefficients.
//' @param withCov Whether any covariates need adjusting.
//' @param type A string indicating the type of outcome, currently either "survival" or "binary".
//' @param seed The seed used for generating the simulation dataset.
//' @return A list of derivatives.
//'
//'@export
//[[Rcpp::export]]
List derivatives(const MatrixXd & Z, const VectorXd & cen,
                 const VectorXd & beta, const VectorXd & y,
                 const MatrixXd & cov, const VectorXd & theta,
                 const int seed, const bool withCov = false,
                 const std::string type = "surv"){

  VectorXd eta;
  VectorXd expeta;

  List r =  {};

  if(withCov)
    eta = Z*beta + cov*theta;
  else
    eta = Z*beta;

  expeta = exp(eta.array());

  if(type == "surv"){
    // calculates eta, u, d, z
    // sorts time in an ascending order
    VectorXd atriskm = revCumSum(expeta).array()*cen.array();
    VectorXd atrisk(atriskm.size());

    for(int i = 0; i < atriskm.size(); i++){
      if(atriskm[i] == 0)
        atrisk[i] = 0;
      else
        atrisk[i] = 1/atriskm[i];
    }

    VectorXd um = expeta.array()*cumSum(atrisk).array();
    VectorXd u = cen - um;

    VectorXd atrisksq = pow(atrisk.array(), 2);



    VectorXd risks = exp(2*eta.array())*cumSum(atrisksq).array();

    VectorXd d = -(u - cen + risks);

    VectorXd z(d.size());

    for(int j = 0; j < eta.size(); j++){
      if(d[j] != 0)
        z[j] = eta[j] + u[j]/d[j];
      else
        z[j] = 0;
    }

    r.insert(0,eta);
    r.insert(1,u);
    r.insert(2,d);
    r.insert(3,z);

  }else{
    // pi must be in the range of (0,1)
    VectorXd nexpeta = exp(- eta.array());
    VectorXd pim = (VectorXd::Ones(eta.size())  + nexpeta);
    VectorXd pi(pim.size());

    for(int j = 0; j < pim.size(); j++){
      pi[j] = 1/pim[j];
    }

    VectorXd dpi = (pi.array() * (1 - pi.array()));
    MatrixXd A = dpi.asDiagonal();

    VectorXd m(eta.size());

    for(int k = 0; k < eta.size(); k++){
      if(A(k,k) != 0)
        m[k] = eta[k] + 1/A(k,k)*(y[k] - pi[k]);
      else
        m[k] = 0;
    }
    r.insert(0,A);
    r.insert(1,m);
  }
  return(r);
}

//' Perform One Iteration of Iterative Re-weighted Least Squares
//'
//' This function performs one iteration of the Iterative Re-weighted Least Squares (IRLS) algorithm based on the given parameters.
//'
//' @param coeff0 Current coefficients for nodes.
//' @param theta0 Current coefficients for covariates that need adjusting.
//' @param cen A vector of censoring statuses or binary outcomes
//' @param y A vector of outcomes.
//' @param Z A matrix of nodes.
//' @param cov A matrix of covariates.
//' @param Cnstrn A matrix of partial ordering constraints.
//' @param coeffLambda A numeric value of the tuning parameter lambda.
//' @param seed The seed used for generating the simulation dataset.
//' @param withCov Whether any covariates need adjusting.
//' @param type A string indicating the type of outcome, either "survival" or "binary".
//'
//' @return A vector of coefficients at the next iteration.
//'
//'@export
//[[Rcpp::export]]
VectorXd IRLSO(const VectorXd & coeff0, const VectorXd & theta0,
               const VectorXd & cen, const VectorXd & y,
               const MatrixXd & Z, const MatrixXd & cov,
               const MatrixXd & Cnstrn, const VectorXd & coeffLambda,
               const int & seed = 0, const bool & withCov = false,
               const string type = "surv"){
  List QP;
  VectorXd V;

  // gets all derivatives
  List dvs = derivatives(Z, cen, coeff0, y, cov, theta0,
                         seed, withCov, type);

  if(type == "surv"){

    VectorXd u0 = dvs.at(1);
    VectorXd d0 = dvs.at(2);
    VectorXd z0 = dvs.at(3);

    if(wInf(u0)||wInf(d0)||wInf(z0)){
      double a = std::numeric_limits<double>::infinity();
      V = VectorXd::Ones(Cnstrn.rows())*a;
      return(V);
    }

    if(withCov){

      MatrixXd d0d = d0.asDiagonal();
      MatrixXd X(Z.rows(), Z.cols()+cov.cols());
      X << Z, cov;
      VectorXd d = (z0.transpose() *d0d * X  - coeffLambda.transpose()).transpose();

      // adds a small value to ensure non-singularity of D
      MatrixXd D = X.transpose()*d0d*X;

      MatrixXd eps = 0.0001*VectorXd::Ones(X.cols()).asDiagonal();
      D = D + eps;

      MatrixXd fullCnstrn(Cnstrn.rows(), Cnstrn.cols()+theta0.size());
      fullCnstrn << Cnstrn, MatrixXd::Zero(Cnstrn.rows(), theta0.size());

      QP = quadprog_solveR(D, d, fullCnstrn.transpose(), VectorXd::Zero(Cnstrn.rows()));

    }else{

      MatrixXd d0d = d0.asDiagonal();
      VectorXd d = (z0.transpose() *d0d * Z  - coeffLambda.transpose()).transpose();

      // adds a small value to ensure non-singularity of D
      MatrixXd D = Z.transpose()*d0d*Z;

      MatrixXd eps = 0.0001*VectorXd::Ones(Z.cols()).asDiagonal();
      D = D + eps;

      QP = quadprog_solveR(D, d, Cnstrn.transpose(), VectorXd::Zero(Cnstrn.rows()));
    }
    V = QP.at(0);
  }else{

    MatrixXd A = dvs.at(0);
    VectorXd m = dvs.at(1);

    if(wInf(m)){
      double a = std::numeric_limits<double>::infinity();
      V = VectorXd::Ones(Cnstrn.rows())*a;
      return(V);
    }

    if(withCov){
      MatrixXd X(Z.rows(), Z.cols()+cov.cols());
      X << Z, cov;
      VectorXd dvec = X.transpose()*A*m - coeffLambda;
      MatrixXd Dmat = X.transpose()*A*X;

      MatrixXd eps = 0.0001*VectorXd::Ones(X.cols()).asDiagonal();
      Dmat = Dmat + eps;

      MatrixXd fullCnstrn(Cnstrn.rows(), Cnstrn.cols()+theta0.size());
      fullCnstrn << Cnstrn, MatrixXd::Zero(Cnstrn.rows(), theta0.size());

      QP = quadprog_solveR(Dmat, dvec, fullCnstrn.transpose(), VectorXd::Zero(Cnstrn.rows()));

    }else{

      VectorXd dvec = Z.transpose()*A*m - coeffLambda;
      MatrixXd Dmat = Z.transpose()*A*Z;


      MatrixXd eps = 0.0001*VectorXd::Ones(Z.cols()).asDiagonal();
      Dmat = Dmat + eps;

      QP = quadprog_solveR(Dmat, dvec, Cnstrn.transpose(), VectorXd::Zero(Cnstrn.rows()));
    }
    V = QP.at(0);
  }
  return(V);
}

//' Iterative Re-weighted Least Squares: All Iterations
//'
//' This function performs all iterations of the Iterative Re-weighted Least Squares (IRLS) algorithm to obtain the Maximum Likelihood Estimates (MLEs) based on the given parameters.
//'
//' @param coeff A vector of initial coefficients for nodes.
//' @param theta A vector of initial coefficients for covariates that need adjusting.
//' @param cen A vector of censoring statuses or binary outcomes.
//' @param y A vector of outcomes.
//' @param Z A matrix of nodes.
//' @param cov A matrix of covariates.
//' @param Cnstrn A matrix of constraints.
//' @param coeffLambda A vector of tuning parameter lambda.
//' @param eps A numeric value denoting the accuracy level.
//' @param maxiter An integer denoting the maximum number of iterations to obtain the MLEs.
//' @param seed The seed used for generating the simulation dataset.
//' @param withCov Whether any covariates need adjusting.
//' @param type A string denoting the outcome type, either "survival" or "binary" outcome.
//' @return The estimates
//'
//'@export
//[[Rcpp::export]]
List IRLSA(const VectorXd & coeff, const VectorXd & theta,
           const VectorXd & cen, const VectorXd & y,
           const MatrixXd & Z, const MatrixXd & cov,
           const MatrixXd & Cnstrn, const VectorXd & coeffLambda,
           const double & eps = 0.0001, const int & maxiter = 10,
           const int & seed = 0, const bool & withCov = false,
           const string type = "surv"){

  double a = std::numeric_limits<double>::infinity();
  double diffCoeff = a;
  double diffTheta;
  List r;

  if(withCov)
    diffTheta = a;

  int nIter = 0;

  VectorXd coff = coeff;
  VectorXd thta = theta;
  VectorXd coeff0;
  VectorXd theta0;


  while(nIter < maxiter && diffCoeff > eps){

    coeff0 = coff;
    theta0 = thta;

    VectorXd beta = IRLSO(coeff0, theta0, cen,
                          y, Z, cov, Cnstrn,
                          coeffLambda, seed,
                          withCov, type);

    if(wInf(beta)){
      if(withCov){
        coff = beta.head(coff.size());
        thta = beta.tail(thta.size());
      }else{
        coff = beta;
      }
      diffCoeff = std::numeric_limits<double>::infinity();
      break;
    }

    if(withCov){

      coff = beta.head(coff.size());

      thta = beta.tail(thta.size());

      diffCoeff = (coff.array() - coeff0.array()).abs().sum()/coff.size();

      diffTheta = (thta.array() - theta0.array()).abs().sum()/thta.size();

      diffCoeff = (diffTheta*thta.size()+diffCoeff*coff.size())/(beta.size());
    }else{
      coff = beta;

      diffCoeff = (coff.array() - coeff0.array()).abs().sum()/coff.size();
    }

    nIter = nIter + 1;
  }

  if(withCov){
    r.insert(0, coff);
    r.insert(1, thta);
    r.insert(2, diffCoeff);
  }else{
    r.insert(0, coff);
    r.insert(1, diffCoeff);
  }
  return(r);
}



//' Perform One Iteration of Iterative Re-weighted Least Squares in the Pruning Step
//'
//' This function performs one iteration of the Iterative Re-weighted Least Squares (IRLS) algorithm in the pruning step based on the given parameters.
//'
//' @param coeff0 Current coefficients for nodes.
//' @param theta0 Current coefficients for covariates that need adjusting.
//' @param cen A vector of censoring statuses or binary outcomes.
//' @param y A vector of outcomes.
//' @param Z A matrix of nodes.
//' @param cov A matrix of covariates.
//' @param LCnstrn A matrix of linear partial ordering constraints.
//' @param qCnstrn A matrix of quadratic constraints.
//' @param mCnstrn An integer denoting how many more constraints are needed.
//' @param minCoeff The lower bound for all the coefficients.
//' @param maxCoeff The upper bound for all the coefficients.
//' @param Lambda A numeric value of the tuning parameter lambda to enforce the quadratic constraints.
//' @param withCov Whether any covariates need adjusting.
//' @param seed The seed used for generating the simulation dataset.
//' @param type A string indicating the type of outcome, either "survival" or "binary" outcome.
//' @return A vector of coefficients at the next iteration.
//'
//'@export
//[[Rcpp::export]]
VectorXd IRLSP(const VectorXd & coeff0, const VectorXd & theta0,
               const VectorXd & cen, const VectorXd & y,
               const MatrixXd & Z, const MatrixXd & cov,
               const MatrixXd & LCnstrn, const MatrixXd & qCnstrn,
               const int & mCnstrn,
               const double & minCoeff, const double & maxCoeff,
               const double & Lambda, const bool & withCov = false,
               const int & seed = 0,
               const string type = "surv"){

  List QP;
  VectorXd V;

  // gets all derivatives
  List dvs = derivatives(Z, cen, coeff0, y, cov, theta0,
                         seed, withCov, type);

  if(type == "surv"){

    VectorXd u0 = dvs.at(1);
    VectorXd d0 = dvs.at(2);
    VectorXd z0 = dvs.at(3);

    if(wInf(u0)||wInf(d0)||wInf(z0)){
      double a = std::numeric_limits<double>::infinity();
      V = VectorXd::Ones(LCnstrn.rows())*a;
      return(V);
    }

    if(withCov){

      MatrixXd d0d = d0.asDiagonal();
      MatrixXd X(Z.rows(), Z.cols()+cov.cols());
      X << Z, cov;
      VectorXd d = (z0.transpose() *d0d * X).transpose();


      MatrixXd D = X.transpose()*d0d*X;
      MatrixXd Q = Lambda*qCnstrn;
      D = D + Q;
      // finds the nearest PD
      D = CnearPD(D);

      MatrixXd fullCnstrn(LCnstrn.rows(), LCnstrn.cols()+theta0.size());
      fullCnstrn << LCnstrn, MatrixXd::Zero(LCnstrn.rows(), theta0.size());

      if(mCnstrn == 0){
        QP = quadprog_solveR(D, d, fullCnstrn.transpose(), VectorXd::Zero(LCnstrn.rows()));
      }else{
        if(mCnstrn == 1){
          VectorXd  bvec(LCnstrn.rows());
          bvec << VectorXd::Zero(LCnstrn.rows()-1), VectorXd::Ones(1)*maxCoeff*(-1);
          QP = quadprog_solveR(D, d, fullCnstrn.transpose(), bvec);
        }else{
          VectorXd  bvec(LCnstrn.rows());
          bvec << VectorXd::Zero(LCnstrn.rows()-2), VectorXd::Ones(1)*maxCoeff*(-1), VectorXd::Ones(1)*minCoeff;
          QP = quadprog_solveR(D, d, fullCnstrn.transpose(), bvec);
        }
      }
    }else{

      MatrixXd d0d = d0.asDiagonal();
      VectorXd d = (z0.transpose() *d0d * Z).transpose();

      MatrixXd D = Z.transpose()*d0d*Z;
      MatrixXd Q = Lambda*qCnstrn;
      D = D + Q;
      // finds the nearest PD
      D = CnearPD(D);

      if(mCnstrn == 0){
        QP = quadprog_solveR(D, d, LCnstrn.transpose(), VectorXd::Zero(LCnstrn.rows()));
      }else{
        if(mCnstrn == 1){
          VectorXd  bvec(LCnstrn.rows());
          bvec << VectorXd::Zero(LCnstrn.rows()-1), VectorXd::Ones(1)*maxCoeff*(-1);
          QP = quadprog_solveR(D, d, LCnstrn.transpose(), bvec);
        }else{
          VectorXd  bvec(LCnstrn.rows());
          bvec << VectorXd::Zero(LCnstrn.rows()-2), VectorXd::Ones(1)*maxCoeff*(-1), VectorXd::Ones(1)*minCoeff;
          QP = quadprog_solveR(D, d, LCnstrn.transpose(), bvec);
        }
      }
    }
    V = QP.at(0);
  }else{
    MatrixXd A = dvs.at(0);
    VectorXd m = dvs.at(1);

    if(wInf(m)){
      double a = std::numeric_limits<double>::infinity();
      V = VectorXd::Ones(LCnstrn.rows())*a;
      return(V);
    }

    if(withCov){
      MatrixXd X(Z.rows(), Z.cols()+cov.cols());
      X << Z, cov;
      VectorXd dvec = X.transpose()*A*m;
      MatrixXd Dmat = X.transpose()*A*X;

      MatrixXd Q = Lambda*qCnstrn;
      Dmat = Dmat + Q;
      // finds the nearest PD
      Dmat = CnearPD(Dmat);

      MatrixXd fullCnstrn(LCnstrn.rows(), LCnstrn.cols()+theta0.size());
      fullCnstrn << LCnstrn, MatrixXd::Zero(LCnstrn.rows(), theta0.size());

      if(mCnstrn == 0){
        QP = quadprog_solveR(Dmat, dvec, fullCnstrn.transpose(), VectorXd::Zero(LCnstrn.rows()));
      }else{
        if(mCnstrn == 1){
          VectorXd  bvec(LCnstrn.rows());
          bvec << VectorXd::Zero(LCnstrn.rows()-1), VectorXd::Ones(1)*maxCoeff*(-1);
          QP = quadprog_solveR(Dmat, dvec, fullCnstrn.transpose(), bvec);
        }else{
          VectorXd  bvec(LCnstrn.rows());
          bvec << VectorXd::Zero(LCnstrn.rows()-2), VectorXd::Ones(1)*maxCoeff*(-1), VectorXd::Ones(1)*minCoeff;
          QP = quadprog_solveR(Dmat, dvec, fullCnstrn.transpose(), bvec);
        }
      }
    }else{

      VectorXd dvec = Z.transpose()*A*m;
      MatrixXd Dmat = Z.transpose()*A*Z;

      MatrixXd Q = Lambda*qCnstrn;
      Dmat = Dmat + Q;
      // finds the nearest PD
      Dmat = CnearPD(Dmat);

      if(mCnstrn == 0){
        QP = quadprog_solveR(Dmat, dvec, LCnstrn.transpose(), VectorXd::Zero(LCnstrn.rows()));
      }else{
        if(mCnstrn == 1){
          VectorXd  bvec(LCnstrn.rows());
          bvec << VectorXd::Zero(LCnstrn.rows()-1), VectorXd::Ones(1)*maxCoeff*(-1);
          QP = quadprog_solveR(Dmat, dvec, LCnstrn.transpose(), bvec);
        }else{
          VectorXd  bvec(LCnstrn.rows());
          bvec << VectorXd::Zero(LCnstrn.rows()-2), VectorXd::Ones(1)*maxCoeff*(-1), VectorXd::Ones(1)*minCoeff;
          QP = quadprog_solveR(Dmat, dvec, LCnstrn.transpose(), bvec);
        }
      }
    }
    V = QP.at(0);
  }

  return(V);
}


//' Perform Parameter Tuning for Each Iteration of Iterative Re-weighted Least Squares in the Pruning Step
//'
//' This function performs parameter tuning for each iteration of the Iterative Re-weighted Least Squares (IRLS) algorithm in the pruning step based on the given parameters.
//'
//' @param coeff A vector of initial coefficients for nodes.
//' @param theta A vector of initial coefficients for covariates that need adjusting.
//' @param cen A vector of censoring statuses or binary outcomes.
//' @param y A vector of outcomes.
//' @param Z A matrix of nodes.
//' @param cov A matrix of covariates.
//' @param LCnstrn A matrix of linear partial ordering constraints.
//' @param qCnstrn A matrix of quadratic constraints.
//' @param mCnstrn An integer denoting how many more constraints are needed.
//' @param minCoeff The lower bound for all the coefficients.
//' @param maxCoeff The upper bound for all the coefficients.
//' @param eps A numeric value denoting the accuracy level.
//' @param maxiter An integer denoting the maximum number of iterations.
//' @param withCov Whether any covariates need adjusting.
//' @param seed The seed used for generating the simulation dataset.
//' @param type A string indicating the type of outcome, either "survival" or "binary" outcome.
//' @return The coefficients after tuning the parameter.
//'
//'@export
// [[Rcpp::export]]

VectorXd IRLSPT(const VectorXd & coeff, const VectorXd & theta,
                const VectorXd & cen, const VectorXd & y,
                const MatrixXd & Z, const MatrixXd & cov,
                const MatrixXd & LCnstrn, const MatrixXd & qCnstrn, const int & mCnstrn,
                const double & minCoeff, const double & maxCoeff,
                const double & eps = 0.0001, const int & maxiter = 10,
                const bool & withCov = false, const int & seed = 0,
                const string type = "surv"){

  double Lambda = 1;
  double a = std::numeric_limits<double>::infinity();
  VectorXd beta_null = VectorXd::Ones(LCnstrn.rows())*a;

  // finds an appropriate lambda when no coefficients go to infinity
  VectorXd beta = IRLSP(coeff, theta,
                        cen, y, Z, cov,
                        LCnstrn, qCnstrn, mCnstrn,
                        minCoeff, maxCoeff,
                        Lambda, withCov,
                        seed, type);

  int k1 = 0;

  while(wInf(beta) && k1 <= 30){

    Lambda = Lambda/2;
    k1 =  k1 + 1;

    beta = IRLSP(coeff, theta,
                 cen, y, Z, cov,
                 LCnstrn, qCnstrn, mCnstrn,
                 minCoeff, maxCoeff,
                 Lambda, withCov,
                 seed, type);
  }

  if(k1 == 31 && wInf(beta)){
    Lambda = 1;
    double k2 = 0;

    while(wInf(beta) && k2 <= 30){

      Lambda = Lambda*2;
      k2 =  k2 + 1;

      beta = IRLSP(coeff, theta,
                   cen, y, Z, cov,
                   LCnstrn, qCnstrn, mCnstrn,
                   minCoeff, maxCoeff,
                   Lambda, withCov,
                   seed, type);
    }

  }

  // no appropriate initial lambda
  if(wInf(beta)){
    return(beta_null);
  }

  // rounds to 3 decimals
  beta = Round(beta,  10*eps);
  double lambda_i = Lambda;

  // finds the maximum and minimum lambda through binary search
  int k3 = 0;
  while((abs(beta.transpose()*qCnstrn*beta) > eps) && (k3 <= 30)){
    Lambda = Lambda*2;
    k3 =  k3 + 1;

    beta = IRLSP(coeff, theta,
                 cen, y, Z, cov,
                 LCnstrn, qCnstrn, mCnstrn,
                 minCoeff, maxCoeff,
                 Lambda, withCov,
                 seed, type);

    if(wInf(beta)){
      break;
    }
    beta = Round(beta,  10*eps);
  }

  double lambda_min, lambda_max;

  if(k3 > 0 && !wInf(beta)){
    if(abs(beta.transpose()*qCnstrn*beta) <= eps){
      lambda_max = Lambda;
      lambda_min = Lambda/2;
    }else{
      return(beta_null);
    }
  }else{
    if(k3 > 0){
      return(beta_null);
    }else{
      int k4 = 0;
      while((abs(beta.transpose()*qCnstrn*beta) <= eps) && (k4 <= 30)){
        Lambda = Lambda/2;
        k4 =  k4 + 1;
        beta = IRLSP(coeff, theta,
                     cen, y, Z, cov,
                     LCnstrn, qCnstrn, mCnstrn,
                     minCoeff, maxCoeff,
                     Lambda, withCov,
                     seed, type);
        if(wInf(beta)){
          Lambda = Lambda*2;
          break;
        }
        beta = Round(beta,  10*eps);
      }
      if(k4 > 0 && !wInf(beta)){
        if(abs(beta.transpose()*qCnstrn*beta) > eps){
          lambda_max = Lambda*2;
          lambda_min = Lambda;
        }else{
          lambda_max = lambda_i;
          lambda_min = Lambda;
        }
      }else{
        if(wInf(beta)){
          lambda_max = lambda_i;
          lambda_min = Lambda;
        }else{
          lambda_max = lambda_i;
          lambda_min = Lambda;
        }
      }
    }
  }

  VectorXd lambda_val;

  if(lambda_min < lambda_max){
    lambda_val = gs(lambda_min, lambda_max);
    beta = IRLSP(coeff, theta,
                 cen, y, Z, cov,
                 LCnstrn, qCnstrn, mCnstrn,
                 minCoeff, maxCoeff,
                 lambda_max, withCov,
                 seed, type);
  }else{
    lambda_val = lambda_min*VectorXd::Ones(1);
    beta = IRLSP(coeff, theta,
                 cen, y, Z, cov,
                 LCnstrn, qCnstrn, mCnstrn,
                 minCoeff, maxCoeff,
                 lambda_min, withCov,
                 seed, type);
  }



  double best_L = logLK(Z, cen, beta.head(Z.cols()), y, cov, beta.tail(cov.cols()),
                        type, withCov);
  double current_L = best_L;

  VectorXd best_beta = beta;

  for (int k = 0; k < lambda_val.size();  k++){

    beta = IRLSP(coeff, theta,
                 cen, y, Z, cov,
                 LCnstrn, qCnstrn, mCnstrn,
                 minCoeff, maxCoeff,
                 lambda_val[k], withCov,
                 seed, type);

    if((!wInf(beta))){
      if((abs(Round(beta,  10*eps).transpose()*qCnstrn*Round(beta,  10*eps)) <= eps)){

        current_L = logLK(Z, cen, beta.head(Z.cols()), y, cov, beta.tail(cov.cols()),
                          type, withCov);

        if(current_L < best_L){
          best_beta = beta;
          best_L = current_L;
        }
      }
    }
  }

  return(best_beta);
}

//' Perform All Iterations of Iterative Re-weighted Least Squares in the Pruning Step
//'
//' This function performs all iterations of the Iterative Re-weighted Least Squares (IRLS) algorithm in the pruning step based on the given parameters.
//'
//' @param coeff A vector of initial coefficients for nodes.
//' @param theta A vector of initial coefficients for covariates that need adjusting.
//' @param cen A vector of censoring statuses or binary outcomes.
//' @param y A vector of outcomes.
//' @param Z A matrix of nodes.
//' @param cov A matrix of covariates.
//' @param LCnstrn A matrix of linear partial ordering constraints.
//' @param qCnstrn A matrix of quadratic constraints.
//' @param mCnstrn An integer denoting how many more constraints are needed.
//' @param eps A numeric value denoting the accuracy level.
//' @param minCoeff The lower bound for all the coefficients.
//' @param maxCoeff The upper bound for all the coefficients.
//' @param maxiter An integer denoting the maximum number of iterations.
//' @param withCov Whether any covariates need adjusting.
//' @param seed The seed used for generating the simulation dataset.
//' @param type A string indicating the type of outcome, either "survival" or "binary" outcome.
//' @return The coefficients under the tuning parameter after several iterations.
//'
//'@export
//[[Rcpp::export]]
List IRLSPAT(const VectorXd & coeff, const VectorXd & theta,
             const VectorXd & cen, const VectorXd & y,
             const MatrixXd & Z, const MatrixXd & cov,
             const MatrixXd & LCnstrn, const MatrixXd & qCnstrn, const int & mCnstrn,
             const double & minCoeff, const double & maxCoeff,
             const double & eps = 0.0001, const int & maxiter = 10,
             const bool & withCov = false, const int & seed = 0,
             const string type = "surv"){

  double a = std::numeric_limits<double>::infinity();

  double diffCoeff = a;
  double diffTheta;
  List r;

  if(withCov)
    diffTheta = a;

  int nIter = 0;

  VectorXd coff = coeff;
  VectorXd thta = theta;
  VectorXd coeff0;
  VectorXd theta0;
  VectorXd beta;



  VectorXd beta_null = VectorXd::Ones(LCnstrn.rows())*a;

  VectorXd best_beta = beta_null;


  double best_L = a;

  double current_L = a;

  while(nIter < maxiter && diffCoeff > eps){

    coeff0 = coff;
    theta0 = thta;

    beta = IRLSPT(coeff0, theta0,
                  cen, y, Z, cov,
                  LCnstrn, qCnstrn, mCnstrn,
                  minCoeff, maxCoeff,
                  eps, maxiter,
                  withCov,
                  seed, type);



    if(wInf(beta)){
      diffCoeff = std::numeric_limits<double>::infinity();
      break;
    }



    if(withCov){

      coff = beta.head(coff.size());

      thta = beta.tail(thta.size());

      diffCoeff = (coff.array() - coeff0.array()).abs().sum()/coff.size();

      diffTheta = (thta.array() - theta0.array()).abs().sum()/thta.size();

      diffCoeff = (diffTheta*thta.size()+diffCoeff*coff.size())/(beta.size());
    }else{
      coff = beta;

      diffCoeff = (coff.array() - coeff0.array()).abs().sum()/coff.size();
    }

    if((!wInf(beta))){

      current_L = logLK(Z, cen, beta.head(Z.cols()), y, cov, beta.tail(cov.cols()),
                        type, withCov);

      if(current_L < best_L){
        best_beta = beta;
        best_L = current_L;
      }

    }

    nIter = nIter + 1;
  }

  r.insert(0, Round(best_beta, 10*eps));
  r.insert(1, diffCoeff);
  r.insert(2, best_L);
  return(r);
}
