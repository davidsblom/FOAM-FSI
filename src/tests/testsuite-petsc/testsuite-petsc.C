
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "gtest/gtest.h"
#include <petscmat.h>
#include <limits.h>

int main(
    int argc,
    char ** argv
    )
{
    ::testing::InitGoogleTest( &argc, argv );
    int result = RUN_ALL_TESTS();
    PetscFinalize();
    return result;
}
