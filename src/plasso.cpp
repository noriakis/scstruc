#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec getEigenValues(arma::mat M) {
    return arma::eig_sym(M);
}

// [[Rcpp::export]]
arma::mat row_replicate (arma::vec weight, int rn) {
	arma::mat rep_mat;
	for (int i=0; i<rn; i++) {
		rep_mat = join_horiz(rep_mat, weight);
	}
	return(rep_mat.t());
}

arma::vec replaceMu (arma::vec weight, double mu) {
	arma::vec new_weights;
	new_weights.zeros(weight.n_elem);
	for (int i = 0; i < weight.n_elem; i++) {
		if (std::fabs(weight[i]) < mu) {
			new_weights[i] <- mu;
		} else {
			new_weights[i] <- weight[i];
		}
	}
	return(new_weights);
}

// [[Rcpp::export]]
arma::mat plasso_fit(arma::mat X, arma::vec y, int maxIter,
	double lambda) {
	
	arma::vec X0, w;
	X0.ones(y.n_elem);
	X = join_horiz(X, X0);
	arma::mat P = X;
	int n = X.n_rows;
	int p = X.n_cols;
	int k_max = 0;
	w.zeros(p);
	arma::mat b = X.t() * y;
	arma::mat G = X.t() * X;
	IntegerVector iter = seq_len(maxIter);
	NumericVector iterNum = as<NumericVector>(iter);
	iterNum = 2 - (8 * (iterNum-1) / maxIter);
	NumericVector muListcpp (iterNum.length());
	for (int i=0; i<iterNum.length(); i++) {
		muListcpp[i] = pow(10, iterNum[i]);
	}
	arma::vec muList = as<arma::vec>(muListcpp);
	
	arma::mat W;
	double mu;
	for (int it=0; it<maxIter; it++) {
		Rcout << "iteration" << it << "\n";
        mu = muList[it];
        mu = lambda;
        w.elem( find_nonfinite(w) ).zeros();
        w = replaceMu(w, mu);
        W = P % row_replicate(w, P.n_rows);
        vec s1 = getEigenValues(W * W.t());
	}
	return 0;
}
