// #include <Rcpp.h> //only include Rcpp.h if RcppArmadillo not included
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double fyxC(double y, arma::mat betmat, arma::colvec X, NumericVector tau) {

  // QR index
  arma::vec xtb = betmat * X;

  arma::vec xtby = xtb - y;

  int ulpos = -99;
  int uupos = -99;

  // find position of one above
  for (int i = 0; i < xtby.size(); i++) {
    if (xtby[i] > 0) {
      uupos = i;
      break;
    }
  }

  // find position of one below
  for (int i = xtby.size(); i > 0; i--) {
    if (xtby[i] <= 0) {
      ulpos = i;
      break;
    }
  }

  //separate code for handling tails
  double lam1 = 1 - tau[0]; 
  double lam2 = tau[tau.size()-1];

  if (ulpos == -99) { //we are way on left tail
    return tau[0] * lam1 * exp(lam1*(-(xtby[0])));
  }

  if (uupos == -99) { // we are way on right tail
    return (1-tau[tau.size()-1]) * lam2 * exp(-lam2 * (-(xtby[tau.size()-1])));
  }

  return (tau[uupos] - tau[ulpos]) / (xtb[uupos] - xtb[ulpos]);
}


// [[Rcpp::export]]
double fvC(double v, int m, NumericVector pi, NumericVector mu, NumericVector sig) {

  double out = 0;

  for (int i; i < m; i++) {
    out += pi[i]/sig[i] * dnorm(NumericVector::create((v-mu[i])/sig[i]))[0]; // * Rcpp::as<double>(dnorm( (v-mu[i])/sig[i] );
  }

  return out;
}

// [[Rcpp::export]]
double fvyxC(double v, arma::mat betmat, int m, NumericVector pi, NumericVector mu, NumericVector sig,
	     double y, arma::colvec x, NumericVector tau) {

  return fyxC(y - v, betmat, x, tau) * fvC(v,m,pi,mu,sig);
}


// [[Rcpp::export]]
NumericVector mh_mcmcC(double startval, int iters, int burnin,
		       double drawsd, arma::mat betmat,
		       int m, NumericVector pi, NumericVector mu,
		       NumericVector sig, double y, arma::mat x,
		       NumericVector tau) {
  //NumericMatrix xx = Rcpp::as<NumericMatrix>(x);
  Environment env = Environment::global_env(); // this will need to be updated
  //Function fv = env["ff.f"];
  Function fvyx = env["fv.yx"];
  NumericVector out(iters);
  out[0] = startval;
  for (int i = 1; i < iters; i++) {
    double trialval = out[i-1] + rnorm(1,0,drawsd)[0];
    double fvold = fvyxC(out[i-1], betmat=betmat, m=m, pi=pi, mu=mu,
			 sig=sig, y=y, x=x, tau=tau);//Rcpp::as<double>(fvyx(out[i-1], betmat=betmat, m=m, pi=pi, mu=mu, sig=sig, y=y, x=x, tau=tau));
    double fvnew = fvyxC(trialval, betmat=betmat, m=m, pi=pi, mu=mu,
			 sig=sig, y=y, x=x, tau=tau);//Rcpp::as<double>(fvyx(trialval, betmat=betmat, m=m, pi=pi, mu=mu, sig=sig, y=y, x=x, tau=tau));
    if (fvnew > fvold) {
      out[i] = trialval;
    } else {
      if ( (fvnew / fvold) >= runif(1)[0]) {
	out[i] = trialval;
      } else {
	out[i] = out[i-1];
      }
    }
  }
  return out;
}
