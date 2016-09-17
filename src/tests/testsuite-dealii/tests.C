
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include <mpi.h>
#include "gtest/gtest.h"
#include <limits.h>
#include <deal.II/base/utilities.h>

int main(int argc,
    char ** argv
    ) {
    dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
