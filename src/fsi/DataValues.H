
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef DataValues_H
#define DataValues_H

#include <Eigen/Dense>

namespace fsi
{
  typedef Eigen::MatrixXd matrix;
  typedef Eigen::VectorXd vector;

  class DataValues
  {
public:

    DataValues();
    ~DataValues();

    void setData( matrix & data );

    void setDataOld( matrix & data );

    matrix data;
    matrix dataprev;
  };
}

#endif
