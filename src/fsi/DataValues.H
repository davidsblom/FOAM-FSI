
/*
 * Copyright [2016] <David Blom>
 */

#ifndef DataValues_H
#define DataValues_H

#include <Eigen/Dense>
#include "fvCFD.H"

namespace fsi {
typedef Eigen::Matrix<scalar, Eigen::Dynamic, Eigen::Dynamic> matrix;
typedef Eigen::Matrix<scalar, Eigen::Dynamic, 1> vector;

class DataValues {
    public:
        DataValues();

        ~DataValues();

        void setData(matrix & data);

        void setDataOld(matrix & data);

        void finalizeTimeStep();

        matrix data;
        matrix dataprev;
        matrix dataPreviousTimeStep;
};
} // namespace fsi

#endif
