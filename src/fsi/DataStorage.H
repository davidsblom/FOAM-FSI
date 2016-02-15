
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#pragma once

#include "SDCSolver.H"
#include "QuadratureInterface.H"

namespace sdc
{
    class DataStorage
    {
public:

        DataStorage(
            std::shared_ptr<fsi::quadrature::IQuadrature<scalar> >,
            int N
            );

        ~DataStorage();

        const fsi::matrix & getFunctions() const;

        const fsi::matrix & getSolutions() const;

        const fsi::matrix integrate(
            const std::vector<scalar> & nodes,
            scalar dt
            ) const;

        const fsi::matrix interpolate(
            const fsi::matrix functions,
            const std::vector<scalar> & nodes
            ) const;

        const fsi::vector getFunction( int substep ) const;

        const fsi::vector getSolution( int substep ) const;

        const fsi::vector getLastSolution() const;

        void initialize(
            int k,
            int N
            );

        void storeFunction(
            const fsi::vector & f,
            int substep
            );

        void storeSolution(
            const fsi::vector & sol,
            int substep
            );

private:

        std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature;
        fsi::matrix F, solStages;
    };
}
