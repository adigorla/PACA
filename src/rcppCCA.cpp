#define ARMA_64BIT_WORD 1

#include <RcppArmadillo.h>
#include <tuple>
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
arma::mat centerCPP(arma::mat& x){
    arma::rowvec xcmeans = mean(x);
    arma::mat xhold(x.n_rows, x.n_cols);
    for(int i = 0 ; i < x.n_rows; i++){
      xhold.row(i) = xcmeans;
    }
    return(x-xhold);
}

// [[Rcpp::export]]
arma::mat getUV(arma::mat& x, arma::mat& eigspc)
{
  arma::mat x_centered = centerCPP(x);
   return(x_centered * eigspc);
 }

// [[Rcpp::export]]
arma::field<arma::mat> doRCCA(arma::mat& x, arma::mat& y){
    Rcout << "Starting doRCCA ..." <<"\n";
	arma::mat Cxy = cov(x, y); // Amat
    //Rcout << "Done with COV XY"<< "\n";
	arma::mat Cxx = cov(x); // Bmat
    //Rcout << "Done with COV X"<< "\n";
	arma::mat Cyy = cov(y); // Cmat
	int p = Cxx.n_rows;
	int q = Cyy.n_rows;
    //Rcout << "The value of p&q : " << p << " "<< q << "\n";
	Cxx = (Cxx + Cxx.t()) / 2.;
	Cyy = (Cyy + Cyy.t()) / 2.;
    //Rcout << "Calculated COV mat"<< "\n";
    //Rcout << "Shape : "<< Cxx.n_rows << " " << Cxx.n_cols << "\n";
    arma::mat Cyy_f;
    arma::mat Cxx_f;
    try {
        Cyy_f = chol(Cyy);

    } catch (std::exception &ex) {
        Rcout << "Cyy not PSD, adding constant to it make PSD"<< "\n";
        Cyy += arma::eye<arma::mat>(Cyy.n_rows,Cyy.n_rows) * 1e-6;
        Cyy_f = chol(Cyy);
    };
    //Rcout << "Done with Cholsky decomp2"<< "\n";
    try {
        Cxx_f = chol(Cxx);

    } catch (std::exception &ex) {
        Rcout << "Cxx not PSD, adding constant to it make PSD"<< "\n";
        Cxx += arma::eye<arma::mat>(Cxx.n_rows,Cxx.n_rows) * 1e-6;
        Cxx_f = chol(Cxx);
    };
    //Rcout << "Done with Cholsky decomp1"<< "\n";
	arma::mat Cxx_fi = arma::inv(Cxx_f);
	arma::mat Cyy_fi = arma::inv(Cyy_f);
    //Rcout << "Done with COV inv"<< "\n";
	arma::mat D = Cxx_fi.t() * Cxy * Cyy_fi;
    //Rcout << "Done with COV matmul"<< "\n";
	arma::mat left;
	arma::mat right;
	arma::vec eigs;
	arma::mat A;
	arma::mat B;
	if (p >= q){
        //Rcout << "in the IF statement"<< "\n";
		svd(left, eigs, right, D);
		A = Cxx_fi * left;
		B = Cyy_fi * right;
	}
	else{
        //Rcout << "in the ELSE statement"<< "\n";
		svd(left, eigs, right, D.t());
		A = Cxx_fi * right;
		B = Cyy_fi * left;
	}
    //Rcout << "Done w/ IF/ELSE!"<< "\n";
	//std::vector<arma::mat> ret_vec;
	//ret_vec.push_back(A);
	//ret_vec.push_back(B);
	arma::field<arma::mat> ret_matrices(2);
	ret_matrices(0) = A;
	ret_matrices(1) = B;
	Rcout << "doRCCA Done!" <<"\n";
	return ret_matrices;
}

// [[Rcpp::export]]
arma::mat getQCPP(arma::mat& x){
  arma::mat Q, R;
  qr_econ(Q, R, x);
  return Q;
}

// [[Rcpp::export]]
arma::field<arma::mat> svdCPP(arma::mat x){
 arma::mat left, right;
 arma::vec eigs;
 svd_econ(left, eigs, right, x);
 arma::field<arma::mat> ret_matrices(3);
 ret_matrices(0) = left;
 ret_matrices(1) = diagmat(eigs);
 ret_matrices(2) = right;
 return ret_matrices;
}

