
/*
 * Copyright [2016] <David Blom>
 */

#ifndef FsiSolver_H
#define FsiSolver_H

#include <deque>
#include <iostream>
#include <list>
#include <memory>
#include "BaseMultiLevelSolver.H"
#include "ConvergenceMeasure.H"
#include "fvCFD.H"

using std::shared_ptr;

namespace fsi {
class FsiSolver {
    public:
        FsiSolver(shared_ptr<BaseMultiLevelSolver> fluid,
            shared_ptr<BaseMultiLevelSolver> solid,
            shared_ptr< std::list<shared_ptr<ConvergenceMeasure> > > convergenceMeasures,
            bool parallel,
            int extrapolationOrder
            );

        virtual ~FsiSolver();

        virtual void evaluate(const vector & input,
            vector & output,
            vector & R
            );

        void extrapolateData();

        void finalizeTimeStep();

        void initTimeStep();

        bool isConvergence();

        bool isConvergence(const vector & x,
            const vector & xprev
            );

        bool isRunning();

        void newMeasurementSeries();

        bool init;
        shared_ptr<BaseMultiLevelSolver> fluid;
        shared_ptr<BaseMultiLevelSolver> solid;
        shared_ptr< std::list<shared_ptr<ConvergenceMeasure> > > convergenceMeasures;
        const bool parallel;
        int N;
        int iter;
        int nbIter;
        bool allConverged;
        vector x;
        const int extrapolationOrder;
        std::deque<vector> previousSolutions;
};
} // namespace fsi

#endif
