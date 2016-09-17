
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "gtest/gtest.h"
#include "SDCTubeFlowFluidSolver.H"
#include "SDCTubeFlowExplicitLinearSolidSolver.H"
#include "GaussRadau.H"
#include "GaussLobatto.H"
#include "Uniform.H"
#include "SDC.H"
#include "SDCFsiSolver.H"
#include "AndersonPostProcessing.H"
#include "RBFCoarsening.H"
#include "RelativeConvergenceMeasure.H"
#include "ResidualRelativeConvergenceMeasure.H"
#include "MinIterationConvergenceMeasure.H"
#include "AitkenPostProcessing.H"
#include "AbsoluteConvergenceMeasure.H"

class SDCFsiExplicitSolidSolverTest : public ::testing::Test {
    protected:
        virtual void SetUp() {
            scalar r0 = 0.2;
            scalar h = 1.0e-3;
            scalar L = 1;
            scalar rho_f = 1.225;
            scalar rho_s = 1.225;
            scalar E0 = 490;
            scalar G = 490;
            scalar nu = 0.5;

            scalar a0 = M_PI * r0 * r0;
            scalar u0 = 0.1;
            scalar p0 = 0;
            scalar dt = 0.001;
            scalar T = 1;
            scalar cmk = std::sqrt(E0 * h / (2 * rho_f * r0));

            int N = 5;
            bool parallel = false;
            int extrapolation = 0;
            scalar tol = 1.0e-3;
            int maxIter = 50;
            scalar initialRelaxation = 1.0e-3;
            int maxUsedIterations = 50;
            int nbReuse = 0;

            scalar singularityLimit = 1.0e-13;
            int reuseInformationStartingFromTimeIndex = 0;
            bool scaling = false;
            bool updateJacobian = false;
            scalar beta = 0.5;
            int minIter = 5;

            std::shared_ptr<tubeflow::SDCTubeFlowFluidSolver> fluid(new tubeflow::SDCTubeFlowFluidSolver(a0, u0, p0, dt, cmk, N, L, T, rho_f));
            std::shared_ptr<fsi::BaseMultiLevelSolver> solid(new tubeflow::SDCTubeFlowExplicitLinearSolidSolver(N, nu, rho_s, h, L, dt, G, E0, r0, p0, T));

            shared_ptr<MultiLevelSolver> fluidSolver(new MultiLevelSolver(fluid, fluid, 0, 0));
            shared_ptr<MultiLevelSolver> solidSolver(new MultiLevelSolver(solid, fluid, 1, 0));

            std::shared_ptr< std::list<std::shared_ptr<ConvergenceMeasure> > > convergenceMeasures;
            convergenceMeasures = std::shared_ptr<std::list<std::shared_ptr<ConvergenceMeasure> > >(new std::list<std::shared_ptr<ConvergenceMeasure> > );

            convergenceMeasures->push_back(std::shared_ptr<ConvergenceMeasure>(new RelativeConvergenceMeasure(0, true, tol)));
            convergenceMeasures->push_back(std::shared_ptr<ConvergenceMeasure>(new MinIterationConvergenceMeasure(0, false, minIter)));

            shared_ptr<MultiLevelFsiSolver> fsi(new MultiLevelFsiSolver(fluidSolver, solidSolver, convergenceMeasures, parallel, extrapolation));

            shared_ptr<PostProcessing> postProcessing(new AndersonPostProcessing(fsi, maxIter, initialRelaxation, maxUsedIterations, nbReuse, singularityLimit, reuseInformationStartingFromTimeIndex, scaling, beta, updateJacobian));

            std::shared_ptr<sdc::SDCFsiSolverInterface> sdcFluidSolver = std::dynamic_pointer_cast<sdc::SDCFsiSolverInterface>(fluid);
            std::shared_ptr<sdc::SDCFsiSolverInterface> sdcSolidSolver = std::dynamic_pointer_cast<sdc::SDCFsiSolverInterface>(solid);

            assert(sdcFluidSolver);
            assert(sdcSolidSolver);

            std::shared_ptr<fsi::SDCFsiSolver> fsiSolver(new fsi::SDCFsiSolver(sdcFluidSolver, sdcSolidSolver, postProcessing, extrapolation));

            int nbNodes = 2;

            std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature;
            quadrature = std::shared_ptr<fsi::quadrature::IQuadrature<scalar> >(new fsi::quadrature::Uniform<scalar>(nbNodes));

            sdc = std::shared_ptr<sdc::SDC> (new sdc::SDC(fsiSolver, quadrature, 1.0e-15, 1, 10));
        }

        virtual void TearDown() {
            sdc.reset();
        }

        std::shared_ptr<sdc::SDC> sdc;
};

TEST_F(SDCFsiExplicitSolidSolverTest, object)
{
    ASSERT_TRUE(true);
}

TEST_F(SDCFsiExplicitSolidSolverTest, timeStep)
{
    sdc->solveTimeStep(0);
    ASSERT_TRUE(sdc->isConverged());
}

TEST_F(SDCFsiExplicitSolidSolverTest, run)
{
    sdc->run();
    ASSERT_TRUE(sdc->isConverged());
}
