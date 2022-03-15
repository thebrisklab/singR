#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// Function that calculates row-wise chordal distance sum
// [[Rcpp::export]]
double chordalD_c(const arma::mat& Ux, const arma::mat& Uy){
  int r = Ux.n_rows;
 
  double dist = r - sum(square(sum(Ux % Uy, 1)));
  
  return dist;
}

// Function that calculates JB
// [[Rcpp::export]]
double calculateJB_c(const arma::mat& U, const arma::mat& X, double alpha = 0.8){
  // Calculate u^{top}X for each u_l and Xj, this is UX
  arma::mat UX = U * X; // r times p

  // Calculate all individual components
  arma::colvec gamma = arma::mean(pow(UX, 3), 1);
  arma::colvec kappa = arma::mean(pow(UX, 4), 1) - 3;
  
  // Calculate JB
  double JB = sum(alpha*square(gamma) + (1-alpha)*square(kappa));
    
  return JB;
}

// Function that calculates joint objective value
// [[Rcpp::export]]
double objectiveJoint_c(const arma::mat& Ux,const arma::mat& Uy,const arma::mat& X, const arma::mat& Y, const arma::mat& normLX, const arma::mat& normLY, double rho, double alpha = 0.8){
  double JBpartX = calculateJB_c(Uy, Y, alpha);
  double JBpartY = calculateJB_c(Ux, X, alpha);
  double Innerproduct = sum(square(sum(normLX % normLY, 1)));
  
  double obj = -JBpartX - JBpartY- 2*rho*Innerproduct;
  return obj;
}

// Function that calculates T(U), JB gradient with respect to U 
// [[Rcpp::export]]
arma::mat calculateT_c(const arma::mat& U, const arma::mat& X, double alpha = 0.8){
  int p = X.n_cols;
  int r = U.n_rows;
  int n = X.n_rows;

  // Calculate u^{top}X for each u_l and Xj, this is UX
  arma::mat UX = U * X;

  // Calculate all individual components
  arma::colvec gamma = arma::mean(pow(UX, 3), 1);
  arma::colvec kappa = arma::mean(pow(UX, 4), 1) - 3;
  arma::mat prod1 = arma::square(UX) * X.t() / p;
  arma::mat prod2 = arma::pow(UX, 3) * X.t() / p;

  // TU must be r by n

  arma::mat TU(r, n);
  TU = 6 * alpha * arma::diagmat(gamma) * prod1 + 8 * (1-alpha) * arma::diagmat(kappa) * prod2;

  //  return(TU)
  return TU;
}

// Funciton that calculates G(U), full gradient with respect to U - one function for X and Y
// [[Rcpp::export]]
arma::mat calculateG_c(const arma::mat&U, const arma::mat& DataW, const arma::mat& invL, const arma::mat& A, double rho, double alpha = 0.8, int r0 = 0){
  
  arma::mat GU = -trans(calculateT_c(U, DataW));
  arma::mat invLU = invL * U.t();
  arma::rowvec normsLU2 = sum(arma::square(invLU), 0);

  if (r0 == 0){
    r0 = U.n_rows;
  }

  double uLa;

  for (int i = 0; i < r0; i++){
    uLa = dot(invLU.col(i), A.row(i));
    GU.col(i) -= 4 * rho * uLa * (invL * trans(A.row(i))/normsLU2[i] - uLa * invL * invLU.col(i)/(std::pow(normsLU2[i],2)));
  }
  return GU;
}


// Full update of both Ux and Uy with chordal distance to measure convergence
// [[Rcpp::export]]
Rcpp::List updateUboth_c(const arma::mat& Ux, const arma::mat& Uy,const arma::mat& xData, const arma::mat& yData,
                      const arma::mat& invLx, const arma::mat& invLy,
                      double rho, int r0, double alpha = 0.8, double tau = 0.01, int maxiter = 1000, double tol = 1e-6){
  // Form standardized LUx and LUy [checked!]
  arma::mat normLX = Ux.rows(0,r0-1) * invLx;
  normLX = diagmat(1/sqrt(sum(arma::square(normLX), 1))) * normLX;
  
  arma::mat normLY = Uy.rows(0,r0-1) * invLy;
  normLY = diagmat(1/sqrt(sum(arma::square(normLY), 1))) * normLY;
  
  // Calculate starting objective value
  double obj1 = objectiveJoint_c(Ux, Uy, xData, yData, normLX, normLY, rho, alpha);
  
  // Initialize error and number of iterations
  int n = Ux.n_cols;
  int rx = Ux.n_rows;
  int ry = Uy.n_rows;
  double error = 100;
  int iter = 0;
  double obj2;
  int tau_ind = 1;
  
  // Store objective values as the algorithm is moving
  arma::colvec fobj(maxiter+1);
  fobj(iter) = obj1;
  
  // Initialize new Ux, Uy; Wx, Wy; GUx, Guy; normLx; normLy
  arma::mat GUx(n, rx);
  arma::mat GUy(n, ry);
  arma::mat Wx(n, n);
  arma::mat Wy(n, n);
  arma::mat Vtaux(rx, n);
  arma::mat Vtauy(ry, n);
  arma::mat normLXtau(r0, n);
  arma::mat normLYtau(r0, n);
  
  // New point starts at supplied starting values
  arma::mat Uxnew = Ux;
  arma::mat Uynew = Uy;
  

  while ((error > tol)&(iter < maxiter)){
    iter++;
    
    // Calculate gradients at current points
    GUx = calculateG_c(Uxnew, xData, invLx, normLY, rho, alpha, r0); // this is n by r
    GUy = calculateG_c(Uynew, yData, invLy, normLX, rho, alpha, r0); // this is n by r

    // Calculate skew-symmetric W
    Wx = GUx * Uxnew - Uxnew.t() * GUx.t();
    Wy = GUy * Uynew - Uynew.t() * GUy.t();

    // Check the value of objective at new Ytau
    obj2 = obj1 + 1;
    tau_ind = 1; // Indicate this is the first time we try that tau

    // Make sure going along descent direction
    while ((obj2 > obj1) & (tau > 1e-14)) {
      // Adjust tau for next search
      if (tau_ind == 0) {tau = 0.8*tau;}

      // Generate new points V(tau) [checked!]
      Vtaux = Uxnew * solve(arma::eye(n,n) + tau*Wx/2, arma::eye(n,n) - tau*Wx/2).t();
      Vtauy = Uynew * solve(arma::eye(n,n) + tau*Wy/2, arma::eye(n,n) - tau*Wy/2).t();

      // Form standardized UX [checked!]
      normLXtau = Vtaux.rows(0,r0-1) * invLx;
      normLXtau = diagmat(1/sqrt(sum(arma::square(normLXtau), 1))) * normLXtau;

      // Form standardized UY [checked!]
      normLYtau = Vtauy.rows(0,r0-1) * invLy;
      normLYtau = diagmat(1/sqrt(sum(arma::square(normLYtau), 1))) * normLYtau;

      // New objective
      obj2 = objectiveJoint_c(Vtaux, Vtauy, xData, yData, normLXtau, normLYtau, rho, alpha);
      
      // Adjust tau indicator if need to return to the loop
      tau_ind = 0;
    }
    
    if (obj2 < obj1){
      obj1 = obj2;
      //error = chordalD_c(Vtaux.rows(0,r0-1), Uxnew.rows(0,r0-1)) + chordalD_c(Vtauy.rows(0,r0-1), Uynew.rows(0,r0-1));
      error = chordalD_c(Vtaux, Uxnew) + chordalD_c(Vtauy, Uynew);
      Uxnew = Vtaux;
      Uynew = Vtauy;
      normLY = normLYtau;
      normLX = normLXtau;
      fobj(iter) = obj1;
    }else{
      error = 0;
    }
  }

  // Otherwise can return the Ux, Uy, normLX, normLY
  return Rcpp::List::create(Rcpp::Named("Ux") = Uxnew, Rcpp::Named("Uy") = Uynew, Rcpp::Named("tau") = tau, Rcpp::Named("iter") = iter, Rcpp::Named("error") = error, Rcpp::Named("fmin") = fobj.subvec(0,iter));
}