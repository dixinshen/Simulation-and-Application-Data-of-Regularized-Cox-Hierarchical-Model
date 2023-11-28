//
// Created by Dixin Shen on 9/26/19.
//

#ifndef PLAYEIGEN_SORT_H
#define PLAYEIGEN_SORT_H

#include <iostream>
#include <Eigen/Dense>
#include <vector>        // std::vector
#include <algorithm>     // std::sort
#include <numeric>       // std:iota

// sort a vector in ascending order, returning ordered indices
std::vector<size_t> sort_index(const Eigen::Ref<const Eigen::VectorXd> &v)
{
    // initialize original index locations
    std::vector<size_t> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);

    // sort index based on comparing values in v
    std::sort(idx.begin(), idx.end(),
         [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

    return idx;
}

// order a matrix based on given indices;
void order(Eigen::MatrixXd &A, std::vector<size_t> &index)
{
    if (A.rows() != index.size()) {
        std::cout << "Error: length of index does not match matrix number of rows!" << std::endl;
    }
    std::vector<Eigen::VectorXd> vec;
    for (int i = 0; i < A.rows(); ++i) {
        vec.emplace_back(A.row(i));
    }

    for (int i = 0; i < A.rows(); ++i)
        A.row(i) = vec[index[i]];

}


#endif //PLAYEIGEN_SORT_H
