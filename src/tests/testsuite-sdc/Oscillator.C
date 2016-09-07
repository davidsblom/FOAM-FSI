
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "Oscillator.H"

Oscillator::Oscillator(int nbTimeSteps,
    scalar dt,
    fsi::vector q0,
    scalar amplitude,
    scalar frequency,
    scalar m,
    scalar k
    )
    :
    sol(q0),
    nbTimeSteps(nbTimeSteps),
    dt(dt),
    q0(q0),
    amplitude(amplitude),
    frequency(frequency),
    m(m),
    k(k),
    t(0) {
    assert(nbTimeSteps > 0);
    assert(dt > 0);
    assert(q0.rows() == 2);
}

void Oscillator::evaluateFunction(const int,
    const fsi::vector & q,
    const scalar t,
    fsi::vector & f
    ) {
    assert(f.rows() == 2);
    assert(q.rows() == 2);

    scalar F = 0.5 * amplitude - 0.5 * amplitude * std::cos(frequency * M_PI * t);

    f(0) = q(1);
    f(1) = F / m - this->k / m * q(0);
}

void Oscillator::finalizeTimeStep() {}

int Oscillator::getDOF() {
    return 2;
}

void Oscillator::getSolution(fsi::vector & solution,
    fsi::vector &
    ) {
    assert(sol.rows() == 2);

    solution = sol;
}

void Oscillator::setSolution(const fsi::vector &,
    const fsi::vector &
    ) {
    assert(false);
}

scalar Oscillator::getEndTime() {
    return nbTimeSteps * dt;
}

scalar Oscillator::getTimeStep() {
    return dt;
}

void Oscillator::nextTimeStep()
{}

void Oscillator::initTimeStep()
{}

void Oscillator::setNumberOfImplicitStages(int)
{}

void Oscillator::implicitSolve(bool,
    const int,
    const int,
    const scalar t,
    const scalar dt,
    const fsi::vector & qold,
    const fsi::vector & rhs,
    fsi::vector & f,
    fsi::vector & result
    ) {
    assert(f.rows() == 2);
    assert(result.rows() == 2);
    assert(qold.rows() == 2);
    assert(rhs.rows() == 2);

    fsi::matrix Ainv(2, 2);
    fsi::vector b(2);
    scalar F;

    Ainv(0, 0) = dt * m / (dt * dt * this->k + m);
    Ainv(0, 1) = dt * dt * m / (dt * dt * this->k + m);
    Ainv(1, 0) = -this->k * dt * dt / (dt * dt * this->k + m);
    Ainv(1, 1) = dt * m / (dt * dt * this->k + m);

    F = 0.5 * amplitude - 0.5 * amplitude * std::cos(frequency * M_PI * t);
    b << qold(0) / dt, F / m + qold(1) / dt;

    b += rhs / dt;
    result = Ainv * b;

    f(0) = result(1);
    f(1) = F / m - this->k / m * result(0);

    sol = result;
    this->t = t;
}
