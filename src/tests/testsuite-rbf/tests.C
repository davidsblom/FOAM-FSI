
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "gtest/gtest.h"
#include <limits.h>
#include "El.hpp"

int main(int argc,
    char ** argv
    ) {
    El::Initialize();
    ::testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    El::Finalize();
    return result;
}
