
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "Piston.H"

Piston::Piston(int nbTimeSteps,
    scalar dt,
    scalar q0,
    scalar qdot0,
    scalar As,
    scalar Ac,
    scalar omega
    )
    :
    SDCSolver(),
    nbTimeSteps(nbTimeSteps),
    dt(dt),
    q0(q0),
    c1(q0 + Ac / std::pow(omega, 2)),
    c2(qdot0 + As / omega),
    As(As),
    Ac(Ac),
    omega(omega),
    N(2),
    q(q0),
    qdot(qdot0),
    t(0),
    timeIndex(0),
    endTime(nbTimeSteps * dt),
    k(0),
    solStages(),
    timeIntegrationScheme(nullptr) {
    assert(nbTimeSteps > 0);
    assert(dt > 0);
}

Piston::Piston(int nbTimeSteps,
    scalar dt,
    scalar q0,
    scalar qdot0,
    scalar As,
    scalar Ac,
    scalar omega,
    std::shared_ptr<sdc::TimeIntegrationScheme> timeIntegrationScheme,
    int k
    )
    :
    SDCSolver(),
    nbTimeSteps(nbTimeSteps),
    dt(dt),
    q0(q0),
    c1(q0 + Ac / std::pow(omega, 2)),
    c2(qdot0 + As / omega),
    As(As),
    Ac(Ac),
    omega(omega),
    N(2),
    q(q0),
    qdot(qdot0),
    t(0),
    timeIndex(0),
    endTime(nbTimeSteps * dt),
    k(k),
    solStages(),
    timeIntegrationScheme(timeIntegrationScheme) {
    assert(nbTimeSteps > 0);
    assert(dt > 0);
    assert(timeIntegrationScheme);
}

Piston::~Piston()
{}

scalar Piston::referenceSolution(scalar t) {
    scalar result = c1;
    result += c2 * t;
    result += -As / std::pow(omega, 2) * std::sin(omega * t);
    result += -Ac / std::pow(omega, 2) * std::cos(omega * t);

    return result;
}

void Piston::evaluateFunction(const int,
    const fsi::vector & q,
    const scalar t,
    fsi::vector & f
    ) {
    assert(f.rows() == 2);

    f(0) = As * std::sin(omega * t);
    f(0) += Ac * std::cos(omega * t);
    f(1) = q(0);
}

void Piston::nextTimeStep() {
    t += dt;
    timeIndex++;

    fsi::vector sol(2);
    sol << qdot, q;

    for (int i = 0; i < k; i++)
        solStages.at(i) = sol;
}

int Piston::getDOF() {
    return 2;
}

scalar Piston::getEndTime() {
    return endTime;
}

void Piston::getSolution(fsi::vector & solution,
    fsi::vector &
    ) {
    assert(solution.rows() == 2);
    solution(0) = qdot;
    solution(1) = q;
}

void Piston::setSolution(const fsi::vector & solution,
    const fsi::vector &
    ) {
    assert(solution.rows() == 2);
    qdot = solution(0);
    q = solution(1);
}

scalar Piston::getTimeStep() {
    return dt;
}

void Piston::run() {
    fsi::vector q(nbTimeSteps + 1), qdot(nbTimeSteps + 1), qold(2), f(2), rhs(2), result(2);

    std::shared_ptr<sdc::SDC> sdc;
    std::shared_ptr<sdc::ESDIRK> esdirk;
    sdc = std::dynamic_pointer_cast<sdc::SDC>(timeIntegrationScheme);
    esdirk = std::dynamic_pointer_cast<sdc::ESDIRK>(timeIntegrationScheme);

    q(0) = -Ac;
    qdot(0) = q(0);

    if (not sdc && not esdirk) {
        setNumberOfImplicitStages(1);

        for (int i = 1; i < nbTimeSteps + 1; i++) {
            nextTimeStep();

            scalar t = dt * i;

            qold << qdot(i - 1), q(i - 1);

            f.setZero();
            rhs.setZero();

            implicitSolve(false, 0, 0, t, dt, qold, rhs, f, result);

            qdot(i) = result(0);
            q(i) = result(1);
        }
    }

    if (sdc || esdirk) {
        setNumberOfImplicitStages(k - 1);

        for (int i = 1; i < nbTimeSteps + 1; i++) {
            nextTimeStep();

            scalar t0 = dt * (i - 1);

            qold << qdot(i - 1), q(i - 1);

            f.setZero();
            rhs.setZero();

            evaluateFunction(0, qold, t0, f);

            timeIntegrationScheme->setOldSolution(i, qold);

            if (sdc)
                sdc->setFunction(-1, f, qold);

            if (esdirk)
                esdirk->setFunction(0, f, qold);

            int nbSweeps = 1;

            if (sdc)
                nbSweeps = 10 * k;

            for (int j = 0; j < nbSweeps; j++) {
                bool corrector = false;

                if (j > 0)
                    corrector = true;

                scalar t = t0;

                int iImplicitStage = 0;

                int nbStages = k - 1;

                if (esdirk && not esdirk->isStageImplicit(esdirk->A(0, 0)))
                    nbStages = k;

                for (int l = 0; l < nbStages; l++) {
                    if (esdirk)
                        if (not esdirk->isStageImplicit(esdirk->A(l, l)))
                            continue;

                    scalar deltaT = 0;

                    if (sdc) {
                        deltaT = dt * sdc->dsdc(l);
                        t += deltaT;
                    }

                    if (esdirk) {
                        deltaT = dt * esdirk->A(l, l);
                        t = t0 + esdirk->C(l) * dt;
                    }

                    Info << "Time = " << t << ", dt = " << dt << ", l = " << l << endl;

                    // run the tests twice to mimic the moving wall velocity boundary
                    // condition.
                    fsi::vector rhs_tmp;

                    for (int i = 0; i < 2; i++) {
                        timeIntegrationScheme->getSourceTerm(corrector, l, j, deltaT, rhs, qold);

                        if (i == 0)
                            rhs_tmp = rhs;

                        if (i == 1) {
                            assert(std::abs(rhs(0) - rhs_tmp(0)) < 1.0e-13);
                            assert(std::abs(rhs(1) - rhs_tmp(1)) < 1.0e-13);
                        }

                        if (sdc)
                            implicitSolve(corrector, l, l, t, deltaT, qold, rhs, f, result);

                        if (esdirk)
                            implicitSolve(corrector, iImplicitStage, 0, t, deltaT, qold, rhs, f, result);

                        timeIntegrationScheme->setFunction(l, f, result);
                    }

                    qdot(i) = result(0);
                    q(i) = result(1);

                    iImplicitStage++;
                }

                timeIntegrationScheme->outputResidual("piston");

                if (timeIntegrationScheme->isConverged() && j >= k)
                    break;
            }
        }
    }
}

void Piston::implicitSolve(bool,
    const int k,
    const int /*kold*/,
    const scalar t,
    const scalar dt,
    const fsi::vector & qold,
    const fsi::vector & rhs,
    fsi::vector & f,
    fsi::vector & result
    ) {
    assert(f.rows() == 2);
    assert(rhs.rows() == 2);
    assert(result.rows() == 2);
    assert(!solStages.empty());

    f(0) = As * std::sin(omega * t);
    f(0) += Ac * std::cos(omega * t);

    // qdot
    result(0) = qold(0) + dt * f(0) + rhs(0);

    // q
    result(1) = qold(1) + std::pow(dt, 2) * f(0) + dt * qold(0) + dt * rhs(0) + rhs(1);

    f(1) = result(0);

    qdot = result(0);
    q = result(1);

    this->t = t;
    this->rhs = rhs;

    solStages.at(k + 1) = result;
}

void Piston::setNumberOfImplicitStages(int k) {
    this->k = k + 1;

    fsi::vector sol(2);
    sol << qdot, q;

    for (int i = 0; i < k + 1; i++)
        solStages.push_back(sol);
}

void Piston::getVariablesInfo(std::deque<int> & dof,
    std::deque<bool> & enabled,
    std::deque<std::string> & names
    ) {
    dof.push_back(1);
    dof.push_back(1);
    enabled.push_back(true);
    enabled.push_back(true);
    names.push_back("qdot");
    names.push_back("q");
}
