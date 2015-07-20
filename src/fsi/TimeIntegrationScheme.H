
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef TimeIntegrationScheme_H
#define TimeIntegrationScheme_H

namespace sdc
{
    class TimeIntegrationScheme
    {
public:

        virtual ~TimeIntegrationScheme(){}

        virtual void getSourceTerm(
            const bool corrector,
            const int k,
            const double deltaT,
            Eigen::VectorXd & rhs,
            Eigen::VectorXd & qold
            ) = 0;

        virtual void setFunction(
            const int k,
            const Eigen::VectorXd & f,
            const Eigen::VectorXd & result
            ) = 0;

        virtual int getNbImplicitStages() = 0;

        virtual void setOldSolution( const Eigen::VectorXd & result ) = 0;

        virtual void outputResidual( std::string name ) = 0;

        virtual bool isConverged() = 0;
    };
}

#endif
