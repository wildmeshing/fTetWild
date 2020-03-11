#pragma once

#include <Eigen/Dense>


namespace floatTetWild {
    void fast_winding_number(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &P, Eigen::VectorXd &W);
}