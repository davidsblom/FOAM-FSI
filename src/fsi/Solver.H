
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef Solver_H
#define Solver_H

#include <ctime>

class Solver
{
public:

    virtual ~Solver(){}

    virtual void run() = 0;

    virtual double timeElapsed()
    {
        time = std::clock() - time;
        return static_cast<double>(time) / CLOCKS_PER_SEC;
    }

    std::clock_t time;
};

#endif
