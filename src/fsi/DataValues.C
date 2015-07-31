
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "DataValues.H"

using namespace fsi;

DataValues::DataValues()
    :
    data(),
    dataprev(),
    dataPreviousTimeStep()
{}

DataValues::~DataValues()
{}

void DataValues::finalizeTimeStep()
{
    dataPreviousTimeStep = data;
}

void DataValues::setData( matrix & data )
{
    this->data = data;

    if ( this->dataprev.cols() == 0 && this->data.cols() > 0 )
    {
        this->dataprev = this->data;
        this->dataprev.setZero();
    }

    if ( dataPreviousTimeStep.cols() == 0 && data.cols() > 0 )
    {
        dataPreviousTimeStep = data;
        dataPreviousTimeStep.setZero();
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

    if ( dataPreviousTimeStep.cols() == 0 && dataprev.cols() > 0 )
    {
        dataPreviousTimeStep = dataprev;
        dataPreviousTimeStep.setZero();
    }
}
