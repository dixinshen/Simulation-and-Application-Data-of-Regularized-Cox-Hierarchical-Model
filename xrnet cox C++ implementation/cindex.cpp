//
// Created by Dixin Shen on 9/25/19.
//

#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>
#include "sort.h"

using namespace Eigen;
using namespace std;
// [[Rcpp::depends(RcppEigen)]]


// [[Rcpp::export]]
double cindex(MatrixXd &y, MatrixXd &x, const Map<VectorXd> &beta)
{
    int N = y.rows();
    vector<size_t> index = sort_index(y.col(0));
    order(y, index);
    order(x, index);
    VectorXd eta = x * beta;
    
    // pairs (i,j) for which the smallest observed time is an event time
    vector<int> wh;
    for (int i = 0; i < N-1; i++) {
        if (y(i,1) == 1) {
            wh.push_back(i);
        }
    }
    double total(0), concordant(0);
    for (auto e : wh) {
        for (int j = e+1; j < N; j++) {
            if (y(j, 0) > y(e, 0)) {
                total += 1;
                if (eta[j] < eta[e]) concordant += 1;
                else if (eta[j] == eta[e]) concordant += 0.5;
            }
        }
    }
    return concordant / total;
}

