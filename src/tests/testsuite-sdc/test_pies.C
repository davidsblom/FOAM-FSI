
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "PIES.H"
#include "Piston.H"
#include "gtest/gtest.h"

class PIESTest : public ::testing::Test {
    protected:
        virtual void SetUp() {
            scalar dt, q0, qdot0, As, Ac, omega, endTime, tol, delta, rho;

            int nbTimeSteps = 100, nbNodes;

            endTime = 100;
            dt = endTime / nbTimeSteps;
            As = 100;
            Ac = As;
            omega = 1;
            q0 = -As;
            qdot0 = -As;
            tol = 1.0e-8;
            delta = 1.0e-9;
            nbNodes = 10;
            rho = 3.15;

            piston = std::shared_ptr<Piston> (new Piston(nbTimeSteps, dt, q0, qdot0, As, Ac, omega));

            pies = std::shared_ptr<sdc::PIES<scalar> > (new sdc::PIES<scalar>(piston, rho, delta, tol, nbNodes, 10 * nbNodes));
        }

        virtual void TearDown() {
            pies.reset();
            piston.reset();
        }

        std::shared_ptr<sdc::PIES<scalar> > pies;
        std::shared_ptr<Piston> piston;
};

TEST_F(PIESTest, coefficients)
{
    ASSERT_TRUE(true);
}

TEST_F(PIESTest, run)
{
    pies->run();

    fsi::vector solution(2), f;
    piston->getSolution(solution, f);

    scalar ref = piston->referenceSolution(100);
    scalar error = std::abs(solution(1) - ref) / std::abs(ref);
    Info << "error = " << error << endl;
    ASSERT_LT(error, 1.0e-8);
}

class PIESRunTest : public ::testing::Test {
    protected:
        virtual void SetUp() {
            scalar dt, q0, qdot0, As, Ac, omega, endTime, tol, delta, rho;

            int nbTimeSteps = 100;

            endTime = 100;
            dt = endTime / nbTimeSteps;
            As = 100;
            Ac = As;
            omega = 1;
            q0 = -As;
            qdot0 = -As;
            tol = 1.0e-8;
            delta = 1.0e-9;
            rho = 3.15;

            pies = std::shared_ptr<sdc::PIES<scalar> > (new sdc::PIES<scalar>(rho, delta, tol));

            piston = std::shared_ptr<Piston> (new Piston(nbTimeSteps, dt, q0, qdot0, As, Ac, omega, pies, pies->nodes.rows()));
        }

        virtual void TearDown() {
            pies.reset();
            piston.reset();
        }

        std::shared_ptr<sdc::PIES<scalar> > pies;
        std::shared_ptr<Piston> piston;
};

TEST_F(PIESRunTest, run)
{
    piston->run();

    fsi::vector solution(2), f;
    piston->getSolution(solution, f);

    scalar ref = piston->referenceSolution(100);
    scalar error = std::abs(solution(1) - ref) / std::abs(ref);
    Info << "error = " << error << endl;
    ASSERT_LT(error, 1.0e-8);
}
