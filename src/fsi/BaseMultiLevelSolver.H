
/*
 * Copyright [2016] <David Blom>
 */

#ifndef BaseMultiLevelSolver_H
#define BaseMultiLevelSolver_H

#include <assert.h>
#include <iostream>
#include <math.h>
#include <memory>
#include "DataValues.H"

namespace fsi
{
    class BaseMultiLevelSolver
    {
        public:
            BaseMultiLevelSolver(
                int N,
                int dim
                )
                :
                N( N ),
                init( false ),
                data( N, dim ),
                initialValue( 0 ),
                substractInitialValue( false ),
                timeIndex( 0 ),
                t( 0 ),
                dim( dim ),
                couplingData()
            {
                assert( N >= 5 );
                assert( dim > 0 );
            }

            BaseMultiLevelSolver(
                int N,
                int dim,
                scalar initialValue
                )
                :
                N( N ),
                init( false ),
                data( N, dim ),
                initialValue( initialValue ),
                substractInitialValue( true ),
                timeIndex( 0 ),
                t( 0 ),
                dim( dim ),
                couplingData()
            {
                assert( N >= 5 );
                assert( dim > 0 );
            }

            virtual ~BaseMultiLevelSolver() {}

            virtual void finalizeTimeStep() = 0;

            virtual void getReadPositions( matrix & readPositions ) = 0;

            virtual void getWritePositions( matrix & writePositions ) = 0;

            virtual void initTimeStep() = 0;

            virtual bool isRunning() = 0;

            virtual void resetSolution() = 0;

            virtual void solve(
                const matrix & input,
                matrix & output
                ) = 0;

            virtual bool interpolateVolField( std::shared_ptr<BaseMultiLevelSolver> )
            {
                return false;
            }

            int N;
            bool init;
            matrix data;
            scalar initialValue;
            bool substractInitialValue;
            int timeIndex;
            scalar t;
            int dim;
            DataValues couplingData;
    };
} // namespace fsi

#endif
