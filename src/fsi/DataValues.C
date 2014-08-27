
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "DataValues.H"

using namespace fsi;

DataValues::DataValues()
{}

DataValues::~DataValues()
{}

void DataValues::setData( matrix & data )
{
  this->data = data;

  if ( this->dataprev.cols() == 0 && this->data.cols() > 0 )
  {
    this->dataprev = this->data;
    this->dataprev.setZero();
  }
}

void DataValues::setDataOld( matrix & data )
{
  this->dataprev = data;

  if ( this->dataprev.cols() == 0 && this->data.cols() > 0 )
  {
    this->dataprev = this->data;
    this->dataprev.setZero();
  }
}
