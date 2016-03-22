
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#ifndef Solver_H
#define Solver_H

#include <ctime>
#include "DataValues.H"

class Solver
{
    public:
        Solver()
            :
            time()
        {}

        virtual ~Solver(){}

        virtual void run() = 0;

        virtual scalar timeElapsed()
        {
            time = std::clock() - time;
            return static_cast<scalar>(time) / CLOCKS_PER_SEC;
        }

        std::clock_t time;
};

#endif
