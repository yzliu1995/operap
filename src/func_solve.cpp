#include<RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
using namespace Rcpp;
using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::MatrixXi;
using Eigen::MatrixXf;
using Eigen::VectorXf;

//' Calculating the sum of a vector
//'
//' This function calculates the sum of a given vector.
//'
//' @param v A numeric vector.
//' @return The sum of the elements in the vector.
//'
//'@export
//[[Rcpp::export]]
double SumL(const VectorXd & v){
  double s = 0;
  for(int i = 0; i < v.size(); i++){
    s+=v[i];
  }
  return(s);
}

//' Calculating the minimum of a vector
//'
//' This function calculates the minimum value of a given vector.
//'
//' @param v A numeric vector.
//' @return The minimum value in the vector.
//'
//'@export
//[[Rcpp::export]]
VectorXd getMin(const VectorXd & v){
  double s = v[0];
  for(int i = 0; i < v.size(); i++){
    if(v[i] < s){
      s = v[i];
    }
  }
  VectorXd mins = VectorXd::Ones(v.size());
  mins = mins*s;
  return(mins);
}


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
List quadprogSolveR(const MatrixXd & Dmat, const VectorXd & dvec,
                    const MatrixXd & Amat)
{
  List result;
  Environment myEnv = Environment::namespace_env("quadprog");
  Function quadR = myEnv["solve.QP"];
  result = quadR(Dmat, dvec, Amat);
  return result;
}

//' Perform Special Summation Over Risk Sets
//'
//' This function calculates a special summation over risk sets based on given parameters.
//'
//' @param indset The index matrix of the risk set.
//' @param value_vec The values to be summed over. In our problem, it represents the exponential of eta.
//' @param ind_vec The vector indicating the failure times to be summed over.
//' @param p The power number.
//'
//' @export
//[[Rcpp::export]]
double sum_risk(const MatrixXd & indset, const VectorXd & value_vec,
                const VectorXd & ind_vec, double p){
  double sum = 0;
  for (int i = 0; i < ind_vec.size(); i++){
    if (ind_vec(i) == 1.) sum += pow(indset.col(i).dot(value_vec), -p);
  }
  return sum;
}

//' Calculate Negative Log Partial Likelihood
//'
//' This function calculates the negative log partial likelihood based on the given parameters.
//'
//' @param X The design matrix.
//' @param indset The risk set matrix.
//' @param beta The vector of coefficients.
//' @param cen The vector of censoring indicators or binary outcomes.
//' @param y The vector of times to failure or binary outcomes.
//' @param type The type of outcome, either "surv" for survival outcome or "bin" for binary outcome.
//'
//' @return The negative log partial likelihood.
//'@export
//[[Rcpp::export]]
double logPL_C(const MatrixXd & X, const MatrixXd & indset,
               const VectorXd & cen, const VectorXd & y,
               const VectorXd & beta,
               const string type = "surv"){
  int n_fail = indset.cols();

  VectorXd eta = X*beta;
  VectorXd expeta = exp(eta.array());
  VectorXd logeta;

  double sum = 0;

  if(type == "bin"){
    logeta = log( (VectorXd::Ones(expeta.size()) + expeta).array() );
    sum  =  -y.transpose() * eta + SumL(logeta);
  }else{
    for (int i = 0; i < n_fail; i++){
      VectorXd col = indset.col(i);
      sum += log(col.dot(expeta));
    }
    sum = sum - cen.dot(eta);
  }

  return(sum);
}

//' Perform Iterative Solution of Group Fused Lasso Penalty
//'
//' This function performs the iterative solution of the group fused lasso penalty based on the given parameters.
//'
//' @param Dmat0 The original quadratic term matrix.
//' @param dvec0 The original linear term vector.
//' @param Apo The partial order matrix.
//' @param lambda The tuning parameter.
//' @param nb The number of boundaries.
//' @param w The weights.
//' @param inibeta The initial beta value for warm start.
//' @param maxiter The maximum number of iterations.
//' @param eps The tolerance number.
//' @param trace Whether to trace each iteration output.
//'
//' @return The solution beta.
//'@export
//[[Rcpp::export]]
VectorXd iter_pen(const MatrixXd & Dmat0, const VectorXd & dvec0,
                  const MatrixXd & Apo,
                  const double lambda, const VectorXd & w,
                  const int nb, const VectorXd & inibeta, const int maxiter,
                  const double eps, const bool trace = false){
  double delta = 1e-8; //small cutoff for group norm in the denominator
  int pqm = Apo.cols();
  if(Apo.rows() != nb || Dmat0.cols() != pqm) Rcout<<"Dimension inconsistent!"<<std::endl;
  VectorXd currbeta = inibeta;
  double diffbeta = 1.;
  int itt = 0;
  VectorXd oldbeta(pqm);
  while(itt < maxiter && diffbeta > eps){
    oldbeta = currbeta;
    //quadratic contribution of the numerator of the penalty
    MatrixXd Q = MatrixXd::Zero(pqm, pqm);
    for(int j = 0; j < nb; j++){
      VectorXd g = VectorXd::Zero(nb);
      g(j) = 1.;
      MatrixXd G = g.replicate(1,1).asDiagonal();
      MatrixXd PGP = Apo.transpose()*G*Apo;
      double norm_g = sqrt(oldbeta.transpose()*PGP*oldbeta);
      //protection from 0 denominator: if norm_g<delta, that PGP would be amplified and force corresponding boundary to 0. Besides, this
      //precision should be higher than the convergence precision.
      if (norm_g > delta)
        Q += (2/norm_g)*PGP*w(j);
      else Q += (2/(delta))*PGP*w(j);
    }

    MatrixXd Dmat = Dmat0 + lambda*Q + 0.01*eps*MatrixXd::Identity(pqm, pqm);
    //re-scales both Dmat and dvec to avoid the overflow of quadprog::solve.QP
    List QP = quadprogSolveR(delta*Dmat, delta*dvec0, Apo.transpose());
    currbeta = QP["solution"];
    diffbeta = (currbeta - oldbeta).norm();

    if(trace){
      Rcpp::Rcout<<"penalty iteration "<<itt<<":"<<std::endl;
      double val = QP["value"];
      Rcpp::Rcout<<"QP value: "<<val<<std::endl;
      Rcpp::Rcout<<"diff beta: "<<diffbeta<<std::endl;
      Rcpp::Rcout<<"current beta: "<<currbeta.transpose()<<std::endl;
      Rcpp::Rcout<<std::endl;
    }

    itt += 1;
  }
  return currbeta;
}

//' Calculate Value of Penalty Term
//'
//' This function calculates the value of the penalty term based on the given parameters.
//'
//' @param Apo The partial order matrix.
//' @param w The adaptive parameter.
//' @param nb The number of boundaries
//' @param beta The vector of coefficients.
//'
//' @return The value of the penalty term.
//'@export
//[[Rcpp::export]]
double val_pen(const MatrixXd Apo, const VectorXd & w,  const int nb, const VectorXd beta){
  if (Apo.rows() != nb) Rcout<<"Inconsistant dimensions of the partial order matrix!"<<std::endl;

  double val = 0;
  for(int j = 0; j<nb; j++){
    VectorXd g = VectorXd::Zero(nb);
    g(j) = 1;
    MatrixXd G = g.replicate(1, 1).asDiagonal();
    double norm = sqrt(beta.transpose()*Apo.transpose()*G*Apo*beta);
    val += w(j)*norm;
  }
  return(val);
}

//' Round Close Values to the Same
//'
//' This function rounds close values in a vector to the same value based on a given tolerance.
//'
//' @param beta The vector of coefficients.
//' @param eps The tolerance for rounding.
//'@export
//[[Rcpp::export]]
VectorXd roundvec(const VectorXd & beta, const double eps){
  int n = beta.size();
  VectorXd beta_round(beta);
  for(int i=0; i<n-1; i++){
    for(int j=i+1; j<n; j++){
      if(std::abs(beta(j) - beta(i))< eps) beta_round(j) = beta_round(i);
    }
  }
  return beta_round;
}


//' Calculate the best coefficients using lasso tree method without parameter tuning
//'
//' This function calculates the best coefficients using lasso tree method without parameter tuning
//'
//' @param X The design matrix.
//' @param indset The risk set.
//' @param Apo The partial ordering matrix.
//' @param nb The number of boundaries.
//' @param cen The vector of censoring indicators or binary outcomes.
//' @param y The vector of times to failure or binary outcomes.
//' @param type The type of outcome, either "surv" for survival outcome or "bin" for binary outcome.
//' @param w The adaptive parameter.
//' @param inibeta The initial beta value for warm start.
//' @param lambda The tuning parameter.
//' @param maxiter The maximum number of iterations.
//' @param eps The tolerance number.
//' @param fullA Boolean variable indicating whether to calculate all elements of the Hessian.
//' @param trace Whether to trace each iteration output.
//'
//' @return A vector of the solution beta.
//'@export
//[[Rcpp::export]]
VectorXd lasso_tree_single(const MatrixXd & X, const MatrixXd & indset,
                           const MatrixXd & Apo,
                           const int nb, const VectorXd & w,
                           const VectorXd & cen,
                           const VectorXd & y,
                           const double lambda, const VectorXd & inibeta,
                           const int maxiter,
                           const double eps, const bool trace,
                           const bool fullA = false,
                           const string type = "surv"){

  int n = X.rows();
  int p = X.cols();
  VectorXd currbeta = inibeta;
  int itt = 0;
  double diffbeta = 1.;
  VectorXd oldbeta(p);
  MatrixXd Dmat(p, p); //quadratic term coefficient in irls
  VectorXd dvec(p); //linear term coefficient in irls
  double oldval, currval, diffval; //tracking the objective values

  currval = logPL_C(X, indset, cen, y, currbeta, type) + lambda*val_pen(Apo, w, nb, currbeta);
  diffval = 1.;


  while(itt < maxiter && diffbeta > eps){
    oldbeta = currbeta;
    oldval = currval;

    VectorXd eta = X*oldbeta;
    VectorXd expeta = exp(eta.array());

    if(type == "surv"){
      VectorXd u(n);
      VectorXd d(n);
      for (int i = 0; i < n; i++){
        u(i) = cen(i) - expeta(i)*sum_risk(indset, expeta, indset.row(i), 1);
        d(i) = u(i) - cen(i) + exp(2*eta(i))*sum_risk(indset, expeta, indset.row(i), 2);
      }
      MatrixXd A(-1.*d.asDiagonal());

      if(fullA){
        for (int i = 0; i < n; i++)
          for (int j = i+1; j < n; j++){
            VectorXd indvec((indset.row(i).array())*(indset.row(j).array()));
            A(j, i) = -2*expeta(i)*expeta(j)*sum_risk(indset, expeta, indvec, 2);
          }
          A = (A + A.transpose())/2.;
      }

      Dmat = X.transpose()*A*X;
      dvec = X.transpose()*(u + A*eta);

      currbeta = iter_pen(Dmat, dvec, Apo, lambda, w, nb, oldbeta, 1, eps, false);
      diffbeta = (currbeta - oldbeta).norm();
    }

    if(type == "bin"){
      // pi must be in the range of (0,1)
      VectorXd nexpeta = exp(- eta.array());
      VectorXd pim = (VectorXd::Ones(eta.size())  + nexpeta);
      VectorXd pi(pim.size());

      for(int j = 0; j < pim.size(); j++){
        pi[j] = 1/pim[j];
      }

      VectorXd dpi = (pi.array() * (1 - pi.array()));
      MatrixXd A = dpi.asDiagonal();

      VectorXd mvec(eta.size());

      for(int k = 0; k < eta.size(); k++){
        if(A(k,k) != 0)
          mvec[k] = eta[k] + 1/A(k,k)*(y[k] - pi[k]);
        else
          mvec[k] = 0;
      }

      Dmat = X.transpose()*A*X;
      dvec = X.transpose()*A*mvec;

      currbeta = iter_pen(Dmat, dvec, Apo, lambda, w, nb, oldbeta, 1, eps, false);
      diffbeta = (currbeta - oldbeta).norm();
    }

    currval = logPL_C(X, indset, cen, y, currbeta, type) + lambda*val_pen(Apo, w, nb, currbeta);
    diffval = currval - oldval;

    if(trace){

      Rcpp::Rcout<<"iteration "<<itt<<":"<<std::endl;

      Rcpp::Rcout<<"diff val:"<<diffval<<std::endl;
      Rcpp::Rcout<<"current objective value: "<<currval<<std::endl;
      Rcpp::Rcout<<"diff beta: "<<diffbeta<<std::endl;
      Rcpp::Rcout<<"current beta: "<<currbeta.transpose()<<std::endl;
      Rcpp::Rcout<<std::endl;
    }

    itt += 1;
  }


  if(itt == maxiter){
    Rcpp::Rcout<<"IRLS doesn't converge! Try increasing maxiter. "<<std::endl;
  }
  currbeta = roundvec(currbeta, eps);
  return currbeta;
}

//' Calculate Best Coefficients using Lasso Tree Method with Parameter Tuning
//'
//' This function calculates the best coefficients using the lasso tree method with parameter tuning.
//'
//' @param X The design matrix.
//' @param indset The risk set.
//' @param Apo The partial ordering matrix.
//' @param nb The number of boundaries.
//' @param cen The vector of censoring indicators or binary outcomes.
//' @param y The vector of times to failure or binary outcomes.
//' @param type The type of outcome, either "surv" for survival outcome or "bin" for binary outcome.
//' @param w The adaptive parameter.
//' @param inibeta The initial beta value for warm start.
//' @param lambda The tuning parameters.
//' @param maxiter The maximum number of iterations.
//' @param eps The tolerance number.
//' @param fullA Boolean variable indicating whether to calculate all elements of the Hessian.
//' @param trace Whether to trace each iteration output.
//'
//' @return A numeric matrix with the first column representing the values of lambda and the remaining columns representing the corresponding coefficients.
//'@export
//[[Rcpp::export]]
MatrixXd lasso_tree_multi(const MatrixXd & X, const MatrixXd & indset,
                          const MatrixXd & Apo,
                          const int nb, const VectorXd & w,
                          const VectorXd & cen, const VectorXd & y,
                          const VectorXd & lambda, const VectorXd & inibeta,
                          const int maxiter, const double eps,
                          const bool trace,
                          const string type = "surv", const bool fullA = false){
  int nlambda = lambda.size();
  int p = X.cols();
  MatrixXd all_beta(nlambda, p);

  VectorXd beta;
  VectorXd beta0 = inibeta;
  for(int i = 0; i < nlambda; i++){
    beta = lasso_tree_single(X, indset, Apo, nb, w, cen, y, lambda(i), beta0, maxiter, eps, trace, fullA, type);
    all_beta.row(i) = beta;
    beta0 = beta;
    if(trace) Rcout<<"lambda = "<<lambda(i)<<", beta's: "<<beta.transpose()<<std::endl;

  }

  return all_beta;
}

