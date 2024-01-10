#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec getEigenValues(arma::mat M) {
    return arma::reverse(arma::eig_sym(M));
}

// [[Rcpp::export]]
arma::mat row_replicate (arma::vec weight, int rn) {
	arma::mat rep_mat;
	for (int i=0; i<rn; i++) {
		rep_mat = join_horiz(rep_mat, weight);
	}
	return(rep_mat.t());
}

// [[Rcpp::export]]
arma::mat col_replicate (arma::vec weight, int rn) {
    arma::mat rep_mat;
    for (int i=0; i<rn; i++) {
        rep_mat = join_horiz(rep_mat, weight);
    }
    return(rep_mat);
}

arma::vec replaceMu (arma::vec weight, double mu) {
	// arma::vec new_weights;
	// new_weights.zeros(weight.n_elem);
	for (unsigned int i = 0; i < weight.n_elem; i++) {
		if (abs(weight[i]) < mu) {
			weight[i] = mu;
		} else {
			weight[i] = weight[i];
		}
	}
	return(weight);
}

// [[Rcpp::export]]
arma::vec plasso_fit(arma::mat X, arma::vec y, int maxIter,
	double lambda, double gamma, double eps) {
	
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
    arma::vec w_old;
	double mu;
	for (int it=0; it<maxIter; it++) {
		Rcout << "lambda " << lambda << ", iteration" << it << "\n";
        mu = muList[it];
        mu = lambda;
        w.elem( find_nonfinite(w) ).zeros();
        w = replaceMu(w, mu);
        Rcout << w << "\n";
        W = P % row_replicate(w, P.n_rows);
        
        vec s1;
        mat U1;
        eig_sym(s1, U1, W*W.t()); 
        U1 = arma::reverse(U1);        
        s1 = replaceMu( arma::reverse(s1), mu );
        s1 = arma::abs(s1);
        Rcout << s1 << "\n";

        W = P % row_replicate(1/w, P.n_rows);
        W.elem( find_nonfinite(W) ).zeros();
        vec s2;
        mat U2;
        eig_sym(s2, U2, W*W.t()); 
        U2 = arma::reverse(U2);  
        s2 = replaceMu( arma::reverse(s2), mu );
        
        s1 = sqrt(s1);
        s2 = sqrt(s2);
        
        U1 = P.t() * U1;
        mat s1Mat = row_replicate( 1 / s1, U1.n_rows );
        s1Mat.elem( find_nonfinite(s1Mat) ).zeros();
        mat D1 = U1 % U1 % s1Mat;
        D1 = arma::sum(D1, 1);
        vec maxi = arma::max(D1);
        vec mini = arma::min(D1);
        D1 = (D1-mini[0]) / (maxi[0]-mini[0]);
        
        U2 = P.t() * U2;
        mat s2Mat = row_replicate( 1 / s2, U2.n_rows );
        s2Mat.elem( find_nonfinite(s2Mat) ).zeros();
        mat D2 = U2 % U2 % s2Mat;
        D2 = arma::sum(D2, 1);
        maxi = arma::max(D2);
        mini = arma::min(D2);
        D2 = (D2-mini[0]) / (maxi[0]-mini[0]);
        
        D1.elem( find_nonfinite(D1) ).zeros();
        D2.elem( find_nonfinite(D2) ).zeros();
        
        mat D = lambda * (gamma * D1 + (1 - gamma) * D2);
        w_old = w;
        
        w = solve(G + diagmat(D), b);
        
        if (sqrt(sum(pow(w - w_old, 2))) < eps) {break;}
	}
	return w;
}
