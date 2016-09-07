
/*
 * Copyright [2016] <David Blom>
 */

#ifndef MultiLevelFsiSolver_H
#define MultiLevelFsiSolver_H

#include <memory>
#include <list>
#include "fvCFD.H"
#include "FsiSolver.H"
#include "MultiLevelSolver.H"

namespace fsi {
class MultiLevelFsiSolver : public FsiSolver {
    public:
        MultiLevelFsiSolver(shared_ptr<MultiLevelSolver> fluidSolver,
            shared_ptr<MultiLevelSolver> solidSolver,
            shared_ptr< std::list<shared_ptr<ConvergenceMeasure> > > convergenceMeasures,
            bool parallel,
            int extrapolationOrder
            );

        ~MultiLevelFsiSolver();

        virtual void evaluate(const fsi::vector & input,
            fsi::vector & output,
            fsi::vector & R
            );

        void initTimeStep();

        bool isConvergence();

        bool isConvergence(const fsi::vector & x,
            const fsi::vector & xprev
            );

        bool isConvergence(shared_ptr< std::list<shared_ptr<ConvergenceMeasure> > > convergenceMeasures);

        void setSurrogateData(fsi::vector & xf,
            matrix & J
            );

        void setUseJacobian(bool useJacobian);

        shared_ptr<MultiLevelSolver> fluidSolver;
        shared_ptr<MultiLevelSolver> solidSolver;

        fsi::vector xf;
        matrix J;
        bool useJacobian;
        int iterCurrentTimeStep;
};
} // namespace fsi

#endif
