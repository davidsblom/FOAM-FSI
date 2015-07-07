
#ifndef QuadratureRules_H
#define QuadratureRules_H

#include <Eigen/Dense>

namespace quadrature {
    void rules(
        std::string rule,
        int n,
        int refine,
        Eigen::VectorXd & nodes,
        Eigen::MatrixXd & smat,
        Eigen::MatrixXd & qmat
        );
}

#endif
