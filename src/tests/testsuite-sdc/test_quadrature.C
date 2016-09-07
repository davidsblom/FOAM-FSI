
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "fvCFD.H"
#include "gtest/gtest.h"
#include "SDC.H"
#include "Piston.H"
#include "Uniform.H"
#include "GaussRadau.H"
#include "GaussLobatto.H"
#include "ClenshawCurtis.H"

using namespace sdc;
using ::testing::TestWithParam;
using ::testing::Bool;
using ::testing::Values;
using ::testing::Combine;

class QuadratureTest : public TestWithParam< std::tr1::tuple<int, std::string> >{
    protected:
        virtual void SetUp() {
            scalar dt, q0, qdot0, As, Ac, omega, endTime, tol;

            int nbNodes = std::tr1::get<0>(GetParam());
            int nbTimeSteps = 100;
            std::string rule = std::tr1::get<1>(GetParam());

            endTime = 100;
            dt = endTime / nbTimeSteps;
            As = 100;
            Ac = As;
            omega = 1;
            q0 = -As;
            qdot0 = -As;
            tol = 1.0e-15;

            std::shared_ptr<sdc::AdaptiveTimeStepper> adaptiveTimeStepper1;
            std::shared_ptr<sdc::AdaptiveTimeStepper> adaptiveTimeStepper2;
            std::shared_ptr<fsi::quadrature::IQuadrature<scalar> > quadrature;

            if (rule == "uniform")
                quadrature = std::shared_ptr<fsi::quadrature::IQuadrature<scalar> >(new fsi::quadrature::Uniform<scalar>(nbNodes));

            if (rule == "gauss-radau")
                quadrature = std::shared_ptr<fsi::quadrature::IQuadrature<scalar> >(new fsi::quadrature::GaussRadau<scalar>(nbNodes));

            if (rule == "clenshaw-curtis")
                quadrature = std::shared_ptr<fsi::quadrature::IQuadrature<scalar> >(new fsi::quadrature::ClenshawCurtis<scalar>(nbNodes));

            if (rule == "gauss-lobatto")
                quadrature = std::shared_ptr<fsi::quadrature::IQuadrature<scalar> >(new fsi::quadrature::GaussLobatto<scalar>(nbNodes));

            adaptiveTimeStepper1 = std::shared_ptr<sdc::AdaptiveTimeStepper> (new sdc::AdaptiveTimeStepper(false));
            adaptiveTimeStepper2 = std::shared_ptr<sdc::AdaptiveTimeStepper> (new sdc::AdaptiveTimeStepper(false));

            piston1 = std::shared_ptr<Piston> (new Piston(nbTimeSteps, dt, q0, qdot0, As, Ac, omega));
            sdc1 = std::shared_ptr<SDC> (new SDC(piston1, quadrature, tol, nbNodes, 10 * nbNodes));

            nbTimeSteps *= 2;
            dt = endTime / nbTimeSteps;

            piston2 = std::shared_ptr<Piston> (new Piston(nbTimeSteps, dt, q0, qdot0, As, Ac, omega));
            sdc2 = std::shared_ptr<SDC> (new SDC(piston2, quadrature, tol, nbNodes, 10 * nbNodes));
        }

        virtual void TearDown() {
            piston1.reset();
            piston2.reset();
            sdc1.reset();
            sdc2.reset();
        }

        std::shared_ptr<Piston> piston1;
        std::shared_ptr<Piston> piston2;
        std::shared_ptr<SDC> sdc1;
        std::shared_ptr<SDC> sdc2;
};

INSTANTIATE_TEST_CASE_P(testParameters, QuadratureTest, ::testing::Combine(Values(2, 4), Values("gauss-radau", "gauss-lobatto", "clenshaw-curtis", "uniform")));

TEST_P(QuadratureTest, order)
{
    sdc1->run();
    sdc2->run();

    fsi::vector solution1(2), solution2(2), f;
    piston1->getSolution(solution1, f);
    piston2->getSolution(solution2, f);

    scalar result1 = solution1(1);
    scalar result2 = solution2(1);

    scalar ref = piston1->referenceSolution(100);
    scalar error1 = std::abs(result1 - ref) / std::abs(ref);
    scalar error2 = std::abs(result2 - ref) / std::abs(ref);

    Info << "error1 = " << error1 << ", error2 = " << error2 << endl;

    int nbTimeSteps = 10;
    scalar order = (std::log10(error1) - std::log10(error2)) / (std::log10(nbTimeSteps * 2) - std::log10(nbTimeSteps));

    int nbNodes = std::tr1::get<0>(GetParam());
    std::string rule = std::tr1::get<1>(GetParam());

    scalar tol = 3.0e-1;

    if (rule == "gauss-radau")
        ASSERT_NEAR(order, nbNodes * 2 - 1, tol);

    if (rule == "gauss-lobatto")
        ASSERT_NEAR(order, nbNodes * 2 - 2, tol);

    if (rule == "clenshaw-curtis")
        ASSERT_NEAR(order, nbNodes, tol);

    if (rule == "uniform")
        ASSERT_NEAR(order, nbNodes, tol);
}
