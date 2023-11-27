//
// Created by Dixin Shen on 9/5/19.
//

#ifndef PLAYEIGEN_QUAD_H
#define PLAYEIGEN_QUAD_H

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>

using namespace Eigen;
using namespace std;

inline void update_quadratic(const Ref<const MatrixXd> &x, const Ref<const VectorXd> &beta, const Ref<const RowVectorXd> &xm,
                             const Ref<const RowVectorXd> &xs, const Ref<const VectorXd> &delta, const Ref<const VectorXi> &ck,
                             const Ref<const VectorXi> &ri, const vector<double> &d, const int &n, const int &m,
                             Ref<VectorXd> W, Ref<VectorXd> r)
{
    VectorXd exp_eta(n);
    for (int i = 0; i < n; i++) {
        exp_eta[i] = exp(x.row(i).cwiseProduct(xs).cwiseProduct(beta.transpose()).sum() - xm.cwiseProduct(xs).cwiseProduct(beta.transpose()).sum());
    }
    double sum_exp_eta_prime = 0;
    VectorXd sum_exp_eta(m);
    for (int i = m-1; i >= 0 && i < m; i--) {
        sum_exp_eta[i] = sum_exp_eta_prime + exp_eta.segment((n-ri[i]), (ri[i]-ri[i+1])).sum();
        sum_exp_eta_prime = sum_exp_eta[i];
    }
    double u_prime = 0;
    double u2_prime = 0;
    for (int k = 0; k < n; k++) {
        if (ck[k+1] == ck[k]) {
            W[k] = exp_eta[k] * u_prime - exp_eta[k] * exp_eta[k] * u2_prime;
            r[k] = (W[k] * log(exp_eta[k]) + delta[k] - exp_eta[k] * u_prime) / n - log(exp_eta[k]) * W[k] / n;
        } else {
            u_prime += d[ck[k+1] - 1] / sum_exp_eta[ck[k+1] - 1];
            u2_prime += d[ck[k+1 - 1]] / (sum_exp_eta[ck[k+1] - 1] * sum_exp_eta[ck[k+1] - 1]);
            W[k] = exp_eta[k] * u_prime - exp_eta[k] * exp_eta[k] * u2_prime;
            r[k] = (W[k] * log(exp_eta[k]) + delta[k] - exp_eta[k] * u_prime) / n - log(exp_eta[k]) * W[k] / n;
        }
    }
}


#endif //PLAYEIGEN_QUAD_H
