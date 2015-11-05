
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <memory>

#include "version.H"
#include "PreciceSolidSolver.H"
#include "SolidSolver.H"

int main(
    int argc,
    char * argv[]
    )
{
    std::shared_ptr<argList> args( new argList( argc, argv ) );

    if ( !args->checkRootCase() )
    {
        FatalError.exit();
    }

    Info << nl << "FOAM-FSI build: " << word( FOAMFSIbuild ) << nl << endl;

    std::shared_ptr<Time> runTime( new Time
        (
            Time::controlDictName,
            args->rootPath(),
            args->caseName()
        ) );

    std::shared_ptr<foamSolidSolver> solid( new SolidSolver( fvMesh::defaultRegion, args, runTime ) );

    PreciceSolidSolver solver( solid );

    solver.run();

    Info << "End\n" << endl;

    label tmp = Pstream::myProcNo();
    reduce( tmp, sumOp<label>() );

    return (0);
}

// ************************************************************************* //
