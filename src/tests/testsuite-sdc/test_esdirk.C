
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <Eigen/Dense>
#include "ESDIRK.H"
#include "Piston.H"
#include "AdaptiveTimeStepper.H"
#include "gtest/gtest.h"
#include "Cos.H"

using namespace sdc;
using ::testing::TestWithParam;
using ::testing::Bool;
using ::testing::Values;
using ::testing::Combine;

/*
 * Implements the piston problem for which a exact solution exists.
 * With the exact solution, the convergence orders can be confirmed
 * exactly for this test case.
 *
 * Reference: A. H. van Zuijlen, Efficient Higher Order Time
 * Integration of Partitioned Systems, PhD thesis, page 22.
 *
 * Notes Gaussian quadrature:
 * http://undergraduate.csse.uwa.edu.au/units/GENG2140/Mech2012/Gaussian_Quadrature.pdf
 * http://en.wikipedia.org/wiki/Gaussian_quadrature
 * http://en.wikipedia.org/wiki/Legendre_polynomials
 */

class ESDIRKTest : public TestWithParam< std::tr1::tuple<int, std::string> >{
    protected:
        virtual void SetUp() {
            scalar dt, q0, qdot0, As, Ac, omega, endTime;

            int nbTimeSteps = std::tr1::get<0>(GetParam());
            std::string method = std::tr1::get<1>(GetParam());

            endTime = 100;
            dt = endTime / nbTimeSteps;
            As = 100;
            Ac = As;
            omega = 1;
            q0 = -As;
            qdot0 = -As;

            std::shared_ptr<sdc::AdaptiveTimeStepper> adaptiveTimeStepper;

            adaptiveTimeStepper = std::shared_ptr<sdc::AdaptiveTimeStepper> (new sdc::AdaptiveTimeStepper(false));

            piston = std::shared_ptr<Piston> (new Piston(nbTimeSteps, dt, q0, qdot0, As, Ac, omega));
            esdirk = std::shared_ptr<ESDIRK> (new ESDIRK(piston, method, adaptiveTimeStepper));

            std::shared_ptr<sdc::ESDIRK> esdirk(new ESDIRK(method));
            piston_esdirk = std::shared_ptr<Piston> (new Piston(nbTimeSteps, dt, q0, qdot0, As, Ac, omega, esdirk, esdirk->getNbImplicitStages() + 1));
        }

        virtual void TearDown() {
            piston.reset();
            esdirk.reset();
            piston_esdirk.reset();
        }

        std::shared_ptr<Piston> piston;
        std::shared_ptr<ESDIRK> esdirk;
        std::shared_ptr<Piston> piston_esdirk;
};

INSTANTIATE_TEST_CASE_P(testParameters, ESDIRKTest, ::testing::Combine(Values(2, 100, 200, 400, 800, 1600, 3200), Values("SDIRK2", "SDIRK2PR", "SDIRK3", "SDIRK4", "ESDIRK3", "ESDIRK4", "ESDIRK5", "ESDIRK53PR", "ESDIRK63PR", "ESDIRK74PR")));

TEST_P(ESDIRKTest, object)
{
    ASSERT_TRUE(true);
}

TEST_P(ESDIRKTest, solveTimeStep)
{
    esdirk->solveTimeStep(0);
}

TEST_P(ESDIRKTest, run)
{
    esdirk->run();

    fsi::vector solution(2), f;
    piston->getSolution(solution, f);

    scalar result = solution(1);
    scalar ref = piston->referenceSolution(100);
    scalar error = std::abs(result - ref) / std::abs(ref);

    int nbTimeSteps = std::tr1::get<0>(GetParam());
    std::string method = std::tr1::get<1>(GetParam());

    ASSERT_NEAR(piston->t, 100, 1.0e-5);

    // ---------------------------------------------

    if (nbTimeSteps == 2 && method == "SDIRK2")
        ASSERT_NEAR(error, 6983.04523052, 1.0e-8);

    if (nbTimeSteps == 100 && method == "SDIRK2")
        ASSERT_NEAR(error, 1.51043910e+00, 1.0e-8);

    if (nbTimeSteps == 200 && method == "SDIRK2")
        ASSERT_NEAR(error, 5.56619516e-01, 1.0e-7);

    if (nbTimeSteps == 400 && method == "SDIRK2")
        ASSERT_NEAR(error, 1.61165601e-01, 1.0e-7);

    if (nbTimeSteps == 800 && method == "SDIRK2")
        ASSERT_NEAR(error, 4.30446278e-02, 1.0e-7);

    if (nbTimeSteps == 1600 && method == "SDIRK2")
        ASSERT_NEAR(error, 1.11061831e-02, 1.0e-7);

    // ---------------------------------------------

    if (nbTimeSteps == 100 && method == "SDIRK3")
        ASSERT_NEAR(error, 2.4490859141990304, 1.0e-8);

    if (nbTimeSteps == 200 && method == "SDIRK3")
        ASSERT_NEAR(error, 0.2897941698557594, 1.0e-7);

    if (nbTimeSteps == 400 && method == "SDIRK3")
        ASSERT_NEAR(error, 0.035416525581987941, 1.0e-7);

    if (nbTimeSteps == 800 && method == "SDIRK3")
        ASSERT_NEAR(error, 0.004382608722869219, 1.0e-7);

    if (nbTimeSteps == 1600 && method == "SDIRK3")
        ASSERT_NEAR(error, 0.00054522799915263077, 1.0e-7);

    // ---------------------------------------------

    if (nbTimeSteps == 100 && method == "ESDIRK3")
        ASSERT_NEAR(error, 2.16401700e-01, 1.0e-8);

    if (nbTimeSteps == 200 && method == "ESDIRK3")
        ASSERT_NEAR(error, 1.59883384e-02, 1.0e-7);

    if (nbTimeSteps == 400 && method == "ESDIRK3")
        ASSERT_NEAR(error, 1.31386639e-03, 1.0e-7);

    if (nbTimeSteps == 800 && method == "ESDIRK3")
        ASSERT_NEAR(error, 1.21257715e-04, 1.0e-7);

    if (nbTimeSteps == 1600 && method == "ESDIRK3")
        ASSERT_NEAR(error, 1.24594213e-05, 1.0e-7);

    // ---------------------------------------------

    if (nbTimeSteps == 100 && method == "ESDIRK4")
        ASSERT_NEAR(error, 1.06552842e-01, 1.0e-8);

    if (nbTimeSteps == 200 && method == "ESDIRK4")
        ASSERT_NEAR(error, 6.62505388e-03, 1.0e-7);

    if (nbTimeSteps == 400 && method == "ESDIRK4")
        ASSERT_NEAR(error, 4.15160887e-04, 1.0e-7);

    if (nbTimeSteps == 800 && method == "ESDIRK4")
        ASSERT_NEAR(error, 2.60153951e-05, 1.0e-7);

    if (nbTimeSteps == 1600 && method == "ESDIRK4")
        ASSERT_NEAR(error, 1.62861015e-06, 1.0e-7);

    // ---------------------------------------------

    if (nbTimeSteps == 100 && method == "ESDIRK5")
        ASSERT_NEAR(error, 9.93536421e-03, 1.0e-8);

    if (nbTimeSteps == 200 && method == "ESDIRK5")
        ASSERT_NEAR(error, 3.04059916e-04, 1.0e-7);

    if (nbTimeSteps == 400 && method == "ESDIRK5")
        ASSERT_NEAR(error, 9.46071494e-06, 1.0e-7);

    if (nbTimeSteps == 800 && method == "ESDIRK5")
        ASSERT_NEAR(error, 2.95451539e-07, 1.0e-7);

    if (nbTimeSteps == 1600 && method == "ESDIRK5")
        ASSERT_NEAR(error, 9.23400561e-09, 1.0e-7);
}

TEST_P(ESDIRKTest, runCompareESDIRK)
{
    piston_esdirk->run();
    esdirk->run();

    fsi::vector solution_piston_esdirk(2), solution_piston(2), f;
    piston->getSolution(solution_piston, f);
    piston_esdirk->getSolution(solution_piston_esdirk, f);

    ASSERT_NEAR(solution_piston_esdirk(0), solution_piston(0), 1.0e-10);
    ASSERT_NEAR(solution_piston_esdirk(1), solution_piston(1), 1.0e-10);
}


TEST(CosTest, ESDIRK)
{
    std::string method = "ESDIRK74PR";

    int nbTimeSteps = 5;
    scalar endTime = 0.05;
    scalar dt = endTime / nbTimeSteps;
    scalar amplitude = 0.2;
    scalar frequency = 5;

    std::shared_ptr<sdc::AdaptiveTimeStepper> adaptiveTimeStepper(new sdc::AdaptiveTimeStepper(false));
    std::shared_ptr<Cos> cos1(new Cos(nbTimeSteps, dt, endTime, amplitude, frequency));
    std::shared_ptr<ESDIRK> esdirk1(new ESDIRK(cos1, method, adaptiveTimeStepper));
    std::shared_ptr<Cos> cos2(new Cos(nbTimeSteps * 2, dt / 2, endTime, amplitude, frequency));
    std::shared_ptr<ESDIRK> esdirk2(new ESDIRK(cos2, method, adaptiveTimeStepper));

    esdirk1->run();
    esdirk2->run();

    scalar ref = 0.5 * amplitude * std::sin(M_PI * frequency * endTime) * M_PI * frequency;
    scalar error1 = std::abs(cos1->f - ref) / std::abs(ref);
    scalar error2 = std::abs(cos2->f - ref) / std::abs(ref);
    scalar order = (std::log10(error1) - std::log10(error2)) / (std::log10(nbTimeSteps * 2) - std::log10(nbTimeSteps));

    std::cout << "error1 = " << error1 << std::endl;
    std::cout << "error2 = " << error2 << std::endl;
    std::cout << "order = " << order << std::endl;
    ASSERT_GE(order, 3.9);

    ref = amplitude * (0.5 - 0.5 * std::cos(M_PI * frequency * endTime));
    error1 = std::abs(cos1->sol - ref) / std::abs(ref);
    error2 = std::abs(cos2->sol - ref) / std::abs(ref);
    order = (std::log10(error1) - std::log10(error2)) / (std::log10(nbTimeSteps * 2) - std::log10(nbTimeSteps));

    std::cout << "error1 = " << error1 << std::endl;
    std::cout << "error2 = " << error2 << std::endl;
    std::cout << "order = " << order << std::endl;
}

class ESDIRKOrderTest : public TestWithParam<std::string>{
    protected:
        virtual void SetUp() {
            scalar dt, q0, qdot0, As, Ac, omega, endTime;

            int nbTimeSteps = 100;
            std::string method = GetParam();

            endTime = 100;
            dt = endTime / nbTimeSteps;
            As = 100;
            Ac = As;
            omega = 1;
            q0 = -As;
            qdot0 = -As;

            std::shared_ptr<sdc::AdaptiveTimeStepper> adaptiveTimeStepper1(new sdc::AdaptiveTimeStepper(false));
            std::shared_ptr<sdc::AdaptiveTimeStepper> adaptiveTimeStepper2(new sdc::AdaptiveTimeStepper(false));

            piston1 = std::shared_ptr<Piston> (new Piston(nbTimeSteps, dt, q0, qdot0, As, Ac, omega));
            esdirk1 = std::shared_ptr<ESDIRK> (new ESDIRK(piston1, method, adaptiveTimeStepper1));

            nbTimeSteps *= 2;
            dt = endTime / nbTimeSteps;

            piston2 = std::shared_ptr<Piston> (new Piston(nbTimeSteps, dt, q0, qdot0, As, Ac, omega));
            esdirk2 = std::shared_ptr<ESDIRK> (new ESDIRK(piston2, method, adaptiveTimeStepper2));
        }

        virtual void TearDown() {
            piston1.reset();
            piston2.reset();
            esdirk1.reset();
            esdirk2.reset();
        }

        std::shared_ptr<Piston> piston1;
        std::shared_ptr<Piston> piston2;
        std::shared_ptr<ESDIRK> esdirk1;
        std::shared_ptr<ESDIRK> esdirk2;
};

INSTANTIATE_TEST_CASE_P(testParameters, ESDIRKOrderTest, ::testing::Values("SDIRK2", "SDIRK2PR", "SDIRK3", "SDIRK4", "ESDIRK3", "ESDIRK4", "ESDIRK5", "ESDIRK53PR", "ESDIRK63PR", "ESDIRK74PR"));

TEST_P(ESDIRKOrderTest, order)
{
    esdirk1->run();
    esdirk2->run();

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

    std::string method = GetParam();

    scalar tol = 0.6;

    if (method == "SDIRK2")
        ASSERT_NEAR(order, 2, tol);

    if (method == "SDIRK2PR")
        ASSERT_NEAR(order, 2, tol);

    if (method == "SDIRK3")
        ASSERT_NEAR(order, 3, tol);

    if (method == "SDIRK4")
        ASSERT_NEAR(order, 4, tol);

    if (method == "ESDIRK3")
        ASSERT_GE(order, 3);

    if (method == "ESDIRK4")
        ASSERT_NEAR(order, 4, tol);

    if (method == "ESDIRK5")
        ASSERT_NEAR(order, 5, tol);

    if (method == "ESDIRK53PR")
        ASSERT_NEAR(order, 3, tol);

    if (method == "ESDIRK63PR")
        ASSERT_NEAR(order, 3, tol);

    if (method == "ESDIRK74PR")
        ASSERT_NEAR(order, 4, tol);
}
