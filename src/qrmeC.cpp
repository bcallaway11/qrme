// #include <Rcpp.h> //only include Rcpp.h if RcppArmadillo not included
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]


//' fyxC
//' Converts quantile regression estimates into density estimates
//'
//' @param y vector of outcomes
//' @param betamat matrix of quantile regression coefficients
//' @param X vector of covariates (should be same dimension as betamat)
//' @param tau vector of values where QR were estimated
//
//' @return value of density at y and X given QR estimates
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

  // some debugging
  // Rcout << "uupos : " << uupos << "\n";
  // Rcout << "ulpos : " << ulpos << "\n";
  // Rcout << "xtb[uupos] : " << xtb[uupos] << "\n";
  // Rcout << "xtb[ulpos] : " << xtb[ulpos] << "\n";
  // Rcout << "tau[uupos] : " << tau[uupos] << "\n";
  // Rcout << "tau[ulpos] : " << tau[ulpos] << "\n";

  //separate code for handling tails
  double lam1 = 1 - tau[0]; 
  double lam2 = tau[tau.size()-1];

  if (ulpos == -99) { //we are way on left tail
    return tau[0] * lam1 * exp(lam1*(-(xtby[0])));
  }

  if (uupos == -99) { // we are way on right tail
    return (1-tau[tau.size()-1]) * lam2 * exp(-lam2 * (-(xtby[tau.size()-1])));
  }

  if ( (tau[uupos] - tau[ulpos]) / (xtb[uupos] - xtb[ulpos]) > 1) {
    return 1; // this is sort of hack, but in a very few number of cases we can get essentially divide by 0s here that
    // result in the same draw over and over...there might be better way to do this though.
  } else {
    return (tau[uupos] - tau[ulpos]) / (xtb[uupos] - xtb[ulpos]);
  }
}

//' fvC
//'
//' Computes the density of the measurement error term given it follows
//'  a mixture of normals distribution
//'
//' @param v value to estimate the density at
//' @param m number of mixture components
//' @param pi vector of mixture probabilities
//' @param mu vector of mixture means
//' @param sig vector of mixture standard deviations
//'
//' @return estimated density of measurement error at v
// [[Rcpp::export]]
double fvC(double v, int m, NumericVector pi, NumericVector mu, NumericVector sig) {

  double out = 0;

  for (int i=0; i < m; i++) {
    out += pi[i]/sig[i] * dnorm(NumericVector::create((v-mu[i])/sig[i]))[0]; 
  }

  return out;
}

//' fvyxC
//'
//' Computes density of measurement error conditional on y and x given
//' QR estimates and distribution of measurement error
//'
//' @param y a particular value of y
//' @param x a vector of x's
//' @inheritParams fvC
//' @inheritParams fyxC
//' @inheritParams mh_mcmcC
//'
//' @return estimate of density of measurement error conditional on y and x
// [[Rcpp::export]]
double fvyxC(double v, arma::mat betmat, int m, NumericVector pi, NumericVector mu, NumericVector sig,
	     double y, arma::colvec x, NumericVector tau) {

  return fyxC(y - v, betmat, x, tau) * fvC(v,m,pi,mu,sig);
}


//' mh_mcmc_innerC
//'
//' Inner part of MCMC algorithm
//'
//' @param y outcome value (for particular observation)
//' @param x vector of covariates (for particular observation)
//' @inheritParams mh_mcmcC
//'
//' @return vector of MCMC draws of measurement error
// [[Rcpp::export]]
NumericVector mh_mcmc_innerC(double startval, int iters, int burnin,
		       double drawsd, arma::mat betmat,
		       int m, NumericVector pi, NumericVector mu,
		       NumericVector sig, double y, arma::mat x,
		       NumericVector tau) {

  NumericVector out(iters);
  out[0] = startval;
  for (int i = 1; i < iters; i++) {
    double trialval = out[i-1] + rnorm(1,0,drawsd)[0];
    double fvold = fvyxC(out[i-1], betmat=betmat, m=m, pi=pi, mu=mu,
    			 sig=sig, y=y, x=x, tau=tau);
    double fvnew = fvyxC(trialval, betmat=betmat, m=m, pi=pi, mu=mu,
    			 sig=sig, y=y, x=x, tau=tau);
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
  return out[Rcpp::Range((burnin), (out.size()-1))];
}


//' imp_sampC
//'
//' return a vector of weights to be used in importance sampling
//' note that, unlike mh_mcmcC, here the measurement error vector
//' has already been drawn and all we need to do is compute weights
//'
//' @param V vector of measurement errors
//' @inheritParams mh_mcmcC
//'
//' @return vector of weights to be used in importance sampling
// [[Rcpp::export]]
NumericVector imp_sampC(NumericVector Y, arma::mat X, NumericVector V, double iters, double drawsd, arma::mat betmat, int m, NumericVector pi, NumericVector mu, NumericVector sig, NumericVector tau) {
  int n = Y.size(); // n includes measurement error so is greater than true number of observations
  double y;
  double v;
  double w1;
  double w2;
  arma::mat x;
  NumericVector weights(n);

  // go ahead and compute denominator of weights
  NumericVector denW = dnorm(V, 0, drawsd);

  // last, compute overall weights as ration of densities
  for (int i = 0; i < n; i++) {
   y = Y[i];
   v = V[i];
   x = X.rows(i,i); // this gets the i-th row
   x = x.t();
   w1 = fvyxC(v, betmat=betmat, m=m, pi=pi, mu=mu,
    			 sig=sig, y=y, x=x, tau=tau);
   w2 = denW[i];
   weights[i] = w1/w2;
  }

  return weights;
}

//' mh_mcmcC
//'
//' @param Y vector of outcomes
//' @param X matrix of covariates
//' @param startval starting value for the markov chain
//' @param iters number of Monte Carlo iterations
//' @param burnin number of first MC iteration to drop
//' @param drawsd the standard deviation for the standard normal draws in the
//'  MH algorithm
//' @param betmat matrix of QR parameters
//' @param m number of mixture components for measurement error
//' @param pi mixture probabilities
//' @param mu means of mixture components
//' @param sig standard deviations of mixture components
//' @param tau which values QR's have been estimated for
// [[Rcpp::export]]
std::vector<double> mh_mcmcC(NumericVector Y, arma::mat X, double startval, int iters,
		   int burnin, double drawsd, arma::mat betmat,
		   int m, NumericVector pi, NumericVector mu,
		   NumericVector sig, NumericVector tau) {

  int n = Y.size();
  int nume = iters-burnin;
  NumericVector e;
  double y;
  arma::mat x;
  std::vector<double> ee;
  ee.reserve(n*nume);

  for (int i = 0; i < n; i++) {
   y = Y[i];
   x = X.rows(i,i); // this gets the i-th row
   x = x.t();
   e = mh_mcmc_innerC(startval=startval, iters=iters, burnin=burnin,
    		      drawsd=drawsd, betmat=betmat,
    		      m=m, pi=pi, mu=mu, sig=sig, y=y, x=x, tau=tau);
   ee.insert(ee.end(), e.begin(), e.end());
  }

  return ee;
  // e = NumericVector::create(y);
  // DataFrame out = DataFrame::create(_["e"]=e);

  // return out;
}

//' fYXmatC
//'
//' Takes n observations of X and a vector of Y's and return
//'   an n x Y.size() matrix of fyx evaluated at those points
//' @inheritParams mh_mcmcC
//' @return matrix of conditional density estimates
// [[Rcpp::export]]
arma::mat fYXmatC(NumericVector Y, arma::mat betmat, arma::mat X, NumericVector tau) {
  int n = X.n_rows;
  int ysize = Y.size();
  arma::mat out(n,ysize);

  for (int i=0; i<n; i++) {
    for (int j=0; j<ysize; j++) {
      out(i,j) = fyxC(Y[j], betmat, X.rows(i,i).t(), tau);
    }
  }

  return out;
}


// [[Rcpp::export]]
class Copula {
  
  char type[20];

//   public:
//   void setType(char t[20]) {
//     type = t;
//   }

// protected:
//   char type[20];

};

//' gumbelCopula
//'
//' Gumbel copula class
// [[Rcpp::export]]
class gumbelCopula: public Copula {
private:
  double alpha;
public:
  gumbelCopula(double a) {
    alpha = a;
  }
  
  std::vector<double> dCopula(std::vector<double> u, std::vector<double> v, double alpha) {

    std::vector<double> out;
    int usize = u.size();
    out.reserve(usize);
    double u1, u2;
    NumericVector otherout(usize);
    
    for (int i=0; i<usize; i++) {
      u1 = u[i];
      u2 = v[i];
      otherout[i] = exp(-pow(pow(-log(u1),alpha) + pow(-log(u2),alpha),(1/alpha))) *  (pow(pow(-log(u1),alpha) + 
										       pow(-log(u2),alpha),((1/alpha) - 1)) * ((1/alpha) * (pow(-log(u2),(alpha - 
										 									  1)) * (alpha * (1/u2)))))
	 * (pow(pow(-log(u1),alpha) + pow(-log(u2),alpha),((1/alpha) - 
																						    					      1)) * ((1/alpha) * (pow(-log(u1),(alpha - 1)) * (alpha * (1/u1)))))
	 - 
	 exp(-pow(pow(-log(u1),alpha) + pow(-log(u2),alpha),(1/alpha))) * (pow(pow(-log(u1),alpha) + 
	 								    pow(-log(u2),alpha),(((1/alpha) - 1) - 1)) * (((1/alpha) - 
	 														   1) * (pow(-log(u2),(alpha - 1)) * (alpha * (1/u2)))) * ((1/alpha) * 
	 																					   (pow(-log(u1),(alpha - 1)) * (alpha * (1/u1)))));
    }

    return(as< std::vector<double> >(otherout));
  }

  double getAlpha() {
    return alpha;
  }
};



//' interpolateC
//'
//' Returns interpolated value at x from parallel arrays ( xData, yData )
//'  Assumes that xData has at least two elements, is sorted and
//' is strictly monotonic increasing
//'
//' @param x vector of x's
//' @param y vector of y's
//' @param xval a particular value to interpolate for
//' @param extrapolate whether or not to linearly extrapolate beyond endpoints
//'  of x
//'
//' @return interpolated value
// [[Rcpp::export]]
double interpolateC(std::vector<double> x, std::vector<double> y, double xval, bool extrapolate) {
   int size = x.size();

   int idx = 0;                                                                  // find left end of interval for interpolation
   if ( xval >= x[size - 2] )                                                 // special case: beyond right end
   {
      idx = size - 2;
   }
   else
   {
      while ( xval > x[idx+1] ) idx++;
   }
   
   double xL = x[idx], yL = y[idx], xR = x[idx+1], yR = y[idx+1];      // points on either side (unless beyond ends)
   if ( !extrapolate )
   {
      if ( xval < xL ) yR = yL;
      if ( xval > xR ) yL = yR;
   }

   // if (xval < xL) yL = ymin;
   // if (xval < xR) yR = ymax;

   

   double dydx = ( yR - yL ) / ( xR - xL );                                    // gradient

   return yL + dydx * ( xval - xL );                                              // linear interpolation
}

//' interpolateMatC
//'
//' vectorized version of interpolateC
//'
//' @param x vector of x's
//' @param ymat matrix of y's
//' @param xval particular value of x to interpolate for
//' @param extrapolate whether or not to extrapolate beyond endpoints of x
//'
//' @return vector of extrapolations
// [[Rcpp::export]]
NumericVector interpolateMatC(std::vector<double> x, arma::mat ymat,
				    double xval, bool extrapolate) {

  int ycols = ymat.n_cols;
  std::vector<double> out;
  out.reserve(ycols);
  std::vector<double> y;
  arma::mat y1;
  NumericMatrix y2;
  NumericVector y3;

  for (int i=0; i<ycols; i++) {
    y1 = ymat.cols(i,i);
    y2 = wrap(y1);
    y3 = y2(_,0);
    y = as< std::vector<double> >(y3);
    out[i] = interpolateC(x, y, xval, extrapolate);
  }

  return wrap(out);
}

// this is not currently used
// leave it here, but not documented
// [[Rcpp::export]]
arma::cube computeFytXC(NumericVector yvals, NumericVector tvals,
			   arma::mat Qyxpreds, arma::mat Ftxpreds,
			   std::vector<double> tau, CharacterVector copula,
			   double copParam) {

  // if (copula=="frank") {
  //   frankCopula cop = frankCopula(copParam);
  // } else if (copula=="gumbel") {
  //   gumbelCopula cop = gumbelCopula(copParam);
  // } else if (copula=="clayton") {
  //   claytonCopula cop = claytonCopula(copParam);
  // } else if (copula=="gaussian") {
  //   normalCopula cop = normalCopula(copParam);
  // }
  gumbelCopula cop = gumbelCopula(copParam);

  int n = Qyxpreds.n_rows;
  int numy = yvals.size();
  int numt = tvals.size();
  NumericVector x3, interpval;
  arma::mat x1, interpval1;
  NumericMatrix x2;
  double y, t, Usum;
  arma::cube out(n,numy,numt);
  NumericVector U = wrap(tau);
  std::vector<double> x, Uinner, copPart;
  x.reserve(Qyxpreds.n_cols);
  Uinner.reserve(U.size());
  copPart.reserve(U.size());
  
  for (int i = 0; i < n; i++) {
    x1 = Qyxpreds.rows(i,i);
    //x1 = x1.t();
    x2 = wrap(x1);
    x3 = x2(0,_);
    x = as< std::vector<double> >(x3);

    for (int j = 0; j < numy; j++) {
      y = yvals[j];

      for (int k = 0; k < numt; k++) {
	t = tvals[k];
	//NumericMatrix Ftx1 = wrap(Ftxpreds.cols(k,k));
	//NumericVector Ftx2 = Ftx1(_,0);
	//std::vector<double> Ftx = as< std::vector<double> >(Ftx2);

	double Ftx1 = Ftxpreds(i,k);
	NumericVector Ftx2(U.size(), Ftx1);
	std::vector<double> Ftx = as< std::vector<double> >(Ftx2);
	
	//interpval = wrap(interpolateMatC(tau, Qyxpreds, 0.5, FALSE));
	//interpval1 = as<arma::mat>(interpval);

	

	Usum = 0;
	copPart = cop.dCopula(tau, Ftx, copParam);

	// some debugging
	// if (i==0 & j==0 & k==0) {
	//   for (int q = 0; q < Ftx.size(); q++) {
	//     Rcpp::Rcout << "Ftx" << std::endl;
	//     Rcpp::Rcout << Ftx[q] << std::endl;
	//     Rcpp::Rcout << "U" << std::endl;
	//     Rcpp::Rcout << tau[q] << std::endl;
	//     Rcpp::Rcout << "copula" << std::endl;
	//     Rcpp::Rcout << copPart[q] << std::endl;
	//     Rcpp::Rcout << "Qyx" << std::endl;
	//     Rcpp::Rcout << x[q] << std::endl;
	//     Rcpp::Rcout << "interpolate" << std::endl;
	//     Rcpp::Rcout << interpolateC(tau, x, U[q], FALSE) << std::endl;  
	//   }
	//   Rcpp::Rcout << "y" << std::endl;
	//   Rcpp::Rcout << y << std::endl;
	// }

	for (int l=0; l<U.size(); l++) {
	  //copPart = cop.dCopula(U, 
	  if ( interpolateC(tau,  x, U[l], FALSE) <= y ) {
	    Usum += copPart[l];
	  }
	}
	out(i,j,k) = Usum/U.size(); //x1*interpval1.t();//x1.t()*interpolateMatC(tau, Qyxmat, 0.5, FALSE); //cop.dCopula(U[1],0.5,copParam); // not sure how to do this...
	}
    }
  }
    //return(interpval);
  return(out);
}

// extra functions for testing copula code
// not used
// [[Rcpp::export]]
NumericVector testCopula(std::vector<double> u, std::vector<double> v, double copParam) {
  gumbelCopula cop = gumbelCopula(copParam);
  NumericVector out(u.size());
  out = wrap(cop.dCopula(u,v,copParam));
  return out;
}
