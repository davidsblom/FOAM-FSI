
#ifndef QuadratureRules_H
#define QuadratureRules_H

#include "BaseMultiLevelSolver.H"

namespace quadrature {
    void rules(
        std::string rule,
        int n,
        int refine,
        fsi::vector & nodes,
        fsi::matrix & smat,
        fsi::matrix & qmat
        );
}

#endif
