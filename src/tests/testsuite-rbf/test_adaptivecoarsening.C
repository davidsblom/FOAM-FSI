
/*
 * Author
 *   David Blom, TU Delft. All rights reserved.
 */

#include "gtest/gtest.h"
#include "ElRBFInterpolation.H"
#include "TPSFunction.H"
#include "RBFFunctionInterface.H"
#include "AdaptiveCoarsening.H"
#include "RBFCoarsening.H"

using namespace rbf;

TEST(AdaptiveCoarsener, interpolate)
{
    std::unique_ptr<RBFFunctionInterface> rbfFunction(new TPSFunction());
    std::unique_ptr<El::DistMatrix<double, El::VR, El::STAR> > positions(new El::DistMatrix<double, El::VR, El::STAR>());
    std::unique_ptr<El::DistMatrix<double, El::VR, El::STAR> > positionsInterpolation(new El::DistMatrix<double, El::VR, El::STAR>());
    std::unique_ptr<El::DistMatrix<double, El::VR, El::STAR> > data(new El::DistMatrix<double, El::VR, El::STAR>());
    El::Zeros(*positionsInterpolation, 8, 1);
    El::Zeros(*positions, 4, 1);
    El::Zeros(*data, positions->Height(), 1);

    const double dx = 1.0 / positions->Height();
    const double dy = 1.0 / double(positionsInterpolation->Height());

    positions->Reserve(positions->LocalHeight() * positions->LocalWidth());

    for (int i = 0; i < positions->LocalHeight(); i++) {
        const int globalRow = positions->GlobalRow(i);

        for (int j = 0; j < positions->LocalWidth(); j++) {
            const int globalCol = positions->GlobalCol(j);
            positions->QueueUpdate(globalRow, globalCol, (globalRow + 1) * dx);
        }
    }

    positionsInterpolation->Reserve(positionsInterpolation->LocalHeight() * positionsInterpolation->LocalWidth());

    for (int i = 0; i < positionsInterpolation->LocalHeight(); i++) {
        const int globalRow = positionsInterpolation->GlobalRow(i);

        for (int j = 0; j < positionsInterpolation->LocalWidth(); j++) {
            const int globalCol = positionsInterpolation->GlobalCol(j);
            positionsInterpolation->QueueUpdate(globalRow, globalCol, (globalRow + 1) * dy);
        }
    }

    AdaptiveCoarsening rbf(1e-3, 1e-2, 2, 10);

    rbf.compute(std::move(rbfFunction), std::move(positions), std::move(positionsInterpolation));

    // Interpolate some data

    data->Reserve(data->LocalHeight() * data->LocalWidth());

    for (int i = 0; i < data->LocalHeight(); i++) {
        const int globalRow = data->GlobalRow(i);

        for (int j = 0; j < data->LocalWidth(); j++) {
            const int globalCol = data->GlobalCol(j);
            const double x = dx * (globalRow + 1);

            data->QueueUpdate(globalRow, globalCol, std::sin(x));
        }
    }

    std::unique_ptr<El::DistMatrix<double, El::VR, El::STAR> > result = rbf.interpolate(std::move(data));
}

TEST(AdaptiveCoarsener, verify_Eigen_coarsening)
{
    std::unique_ptr<RBFFunctionInterface> rbfFunction(new TPSFunction());
    std::unique_ptr<El::DistMatrix<double, El::VR, El::STAR> > positions(new El::DistMatrix<double, El::VR, El::STAR>());
    std::unique_ptr<El::DistMatrix<double, El::VR, El::STAR> > positionsInterpolation(new El::DistMatrix<double, El::VR, El::STAR>());
    std::unique_ptr<El::DistMatrix<double, El::VR, El::STAR> > data(new El::DistMatrix<double, El::VR, El::STAR>());
    El::Zeros(*positionsInterpolation, 20, 1);
    El::Zeros(*positions, 10, 1);
    El::Zeros(*data, positions->Height(), 1);
    matrix positionsEigen(positions->Height(), positions->Width()), positionsInterpolationEigen(positionsInterpolation->Height(), positionsInterpolation->Width());
    std::shared_ptr<RBFFunctionInterface> rbfFunctionEigen(new TPSFunction());
    matrix valuesEigen(positionsEigen.rows(), positionsEigen.cols()), valuesInterpolationEigen;

    const double dx = 1.0 / positions->Height();
    const double dy = 1.0 / double(positionsInterpolation->Height());

    for (int i = 0; i < positionsEigen.rows(); i++)
        positionsEigen(i, 0) = (i + 1) * dx;

    for (int i = 0; i < positionsInterpolationEigen.rows(); i++)
        positionsInterpolationEigen(i, 0) = (i + 1) * dy;

    for (int i = 0; i < valuesEigen.rows(); i++)
        valuesEigen(i, 0) = std::sin((i + 1) * dx);

    std::shared_ptr<RBFInterpolation> rbfInterpolator(new RBFInterpolation(rbfFunctionEigen, false, true));
    RBFCoarsening rbfEigen(rbfInterpolator, true, true, false, 1e-3, 0.1, 2, 5, false);

    rbfEigen.compute(positionsEigen, positionsInterpolationEigen);
    rbfEigen.interpolate(valuesEigen, valuesInterpolationEigen);

    positions->Reserve(positions->LocalHeight() * positions->LocalWidth());

    for (int i = 0; i < positions->LocalHeight(); i++) {
        const int globalRow = positions->GlobalRow(i);

        for (int j = 0; j < positions->LocalWidth(); j++) {
            const int globalCol = positions->GlobalCol(j);
            positions->QueueUpdate(globalRow, globalCol, (globalRow + 1) * dx);
        }
    }

    positionsInterpolation->Reserve(positionsInterpolation->LocalHeight() * positionsInterpolation->LocalWidth());

    for (int i = 0; i < positionsInterpolation->LocalHeight(); i++) {
        const int globalRow = positionsInterpolation->GlobalRow(i);

        for (int j = 0; j < positionsInterpolation->LocalWidth(); j++) {
            const int globalCol = positionsInterpolation->GlobalCol(j);
            positionsInterpolation->QueueUpdate(globalRow, globalCol, (globalRow + 1) * dy);
        }
    }

    positions->ProcessQueues();
    positionsInterpolation->ProcessQueues();

    std::vector<double> bufferPos;
    positions->ReservePulls(positions->Height() * positions->Width());

    for (int i = 0; i < positions->Height(); i++) {
        for (int j = 0; j < positions->Width(); j++) {
            positions->QueuePull(i, j);
        }
    }

    positions->ProcessPullQueue(bufferPos);

    int index = 0;

    for (int i = 0; i < positions->Height(); i++) {
        for (int j = 0; j < positions->Width(); j++) {
            EXPECT_NEAR(positionsEigen(i, j), bufferPos[index], 1e-13);
            index++;
        }
    }

    positionsInterpolation->ReservePulls(positionsInterpolation->Height() * positionsInterpolation->Width());
    bufferPos.clear();

    for (int i = 0; i < positionsInterpolation->Height(); i++) {
        for (int j = 0; j < positionsInterpolation->Width(); j++) {
            positionsInterpolation->QueuePull(i, j);
        }
    }

    positionsInterpolation->ProcessPullQueue(bufferPos);

    index = 0;

    for (int i = 0; i < positionsInterpolation->Height(); i++) {
        for (int j = 0; j < positionsInterpolation->Width(); j++) {
            EXPECT_NEAR(positionsInterpolationEigen(i, j), bufferPos[index], 1e-13);
            index++;
        }
    }

    AdaptiveCoarsening rbf(1e-3, 1e-2, 2, 5);
    rbf.compute(std::move(rbfFunction), std::move(positions), std::move(positionsInterpolation));

    // Interpolate some data

    data->Reserve(data->LocalHeight() * data->LocalWidth());

    for (int i = 0; i < data->LocalHeight(); i++) {
        const int globalRow = data->GlobalRow(i);

        for (int j = 0; j < data->LocalWidth(); j++) {
            const int globalCol = data->GlobalCol(j);
            const double x = dx * (globalRow + 1);

            data->QueueUpdate(globalRow, globalCol, std::sin(x));
        }
    }

    bufferPos.clear();
    data->ReservePulls(data->Height() * data->Width());

    for (int i = 0; i < data->Height(); i++) {
        for (int j = 0; j < data->Width(); j++) {
            data->QueuePull(i, j);
        }
    }

    data->ProcessPullQueue(bufferPos);

    index = 0;

    for (int i = 0; i < data->Height(); i++) {
        for (int j = 0; j < data->Width(); j++) {
            EXPECT_NEAR(valuesEigen(i, j), bufferPos[index], 1e-13);
            index++;
        }
    }

    std::unique_ptr<El::DistMatrix<double, El::VR, El::STAR> > result = rbf.interpolate(std::move(data));

    std::vector<double> buffer;
    result->ReservePulls(result->Height() * result->Width());

    for (int i = 0; i < result->Height(); i++)
        for (int j = 0; j < result->Width(); j++)
            result->QueuePull(i, j);

    result->ProcessPullQueue(buffer);

    index = 0;

    for (int i = 0; i < result->Height(); i++) {
        for (int j = 0; j < result->Width(); j++) {
            ASSERT_NEAR(valuesInterpolationEigen(i, j), buffer[index], 1e-13);
            index++;
        }
    }
}

TEST(AdaptiveCoarsener, 2d_grid)
{
    std::unique_ptr<RBFFunctionInterface> rbfFunction(new TPSFunction());
    std::unique_ptr<El::DistMatrix<double, El::VR, El::STAR> > positions(new El::DistMatrix<double, El::VR, El::STAR>());
    std::unique_ptr<El::DistMatrix<double, El::VR, El::STAR> > positionsInterpolation(new El::DistMatrix<double, El::VR, El::STAR>());
    std::unique_ptr<El::DistMatrix<double, El::VR, El::STAR> > data(new El::DistMatrix<double, El::VR, El::STAR>());
    El::Zeros(*positions, 20, 2);
    El::Zeros(*positionsInterpolation, 50, 2);
    El::Zeros(*data, positions->Height(), 1);

    matrix positionsEigen(positions->Height(), positions->Width()), positionsInterpolationEigen(positionsInterpolation->Height(), positionsInterpolation->Width());
    std::shared_ptr<RBFFunctionInterface> rbfFunctionEigen(new TPSFunction());
    matrix valuesEigen(positionsEigen.rows(), 1), valuesInterpolationEigen;

    positions->Reserve(positions->LocalHeight() * positions->LocalWidth());

    for (int i = 0; i < positions->LocalHeight(); i++) {
        const int globalRow = positions->GlobalRow(i);

        for (int j = 0; j < positions->LocalWidth(); j++) {
            const int globalCol = positions->GlobalCol(j);
            positions->QueueUpdate(globalRow, globalCol, ((double) rand() / (RAND_MAX)) * 4.0 - 2.0);
        }
    }

    positionsInterpolation->Reserve(positionsInterpolation->LocalHeight() * positionsInterpolation->LocalWidth());

    for (int i = 0; i < positionsInterpolation->LocalHeight(); i++) {
        const int globalRow = positionsInterpolation->GlobalRow(i);

        for (int j = 0; j < positionsInterpolation->LocalWidth(); j++) {
            const int globalCol = positionsInterpolation->GlobalCol(j);
            int row = globalRow % 10;
            int col = globalRow / 10;

            if (globalCol == 0)   // x-coordinate
                positionsInterpolation->QueueUpdate(globalRow, globalCol, -2.0 + 0.05 * row);
            else // y-coordinate
                positionsInterpolation->QueueUpdate(globalRow, globalCol, -2.0 + 0.05 * col);
        }
    }

    positions->ProcessQueues();
    positionsInterpolation->ProcessQueues();

    std::vector<double> bufferPos;
    positions->ReservePulls(positions->Height() * positions->Width());

    for (int i = 0; i < positions->Height(); i++) {
        for (int j = 0; j < positions->Width(); j++) {
            positions->QueuePull(i, j);
        }
    }

    positions->ProcessPullQueue(bufferPos);

    int index = 0;

    for (int i = 0; i < positions->Height(); i++) {
        for (int j = 0; j < positions->Width(); j++) {
            positionsEigen(i, j) = bufferPos[index];
            index++;
        }
    }

    // data

    data->Reserve(data->LocalHeight() * data->LocalWidth());

    for (int i = 0; i < data->LocalHeight(); i++) {
        const int globalRow = data->GlobalRow(i);

        const double x = bufferPos[globalRow * 2];
        const double y = bufferPos[globalRow * 2 + 1];

        data->QueueUpdate(globalRow, 0, x * std::exp(-x * x - y * y));
    }

    data->ProcessQueues();

    positionsInterpolation->ReservePulls(positionsInterpolation->Height() * positionsInterpolation->Width());
    bufferPos.clear();

    for (int i = 0; i < positionsInterpolation->Height(); i++) {
        for (int j = 0; j < positionsInterpolation->Width(); j++) {
            positionsInterpolation->QueuePull(i, j);
        }
    }

    positionsInterpolation->ProcessPullQueue(bufferPos);

    index = 0;

    for (int i = 0; i < positionsInterpolation->Height(); i++) {
        for (int j = 0; j < positionsInterpolation->Width(); j++) {
            positionsInterpolationEigen(i, j) = bufferPos[index];
            index++;
        }
    }

    AdaptiveCoarsening rbf(1e-3, 1e-2, 2, 15);

    rbf.compute(std::move(rbfFunction), std::move(positions), std::move(positionsInterpolation));

    std::unique_ptr<El::DistMatrix<double, El::VR, El::STAR> > result = rbf.interpolate(data);

    bufferPos.clear();
    data->ReservePulls(data->Height() * data->Width());

    for (int i = 0; i < data->Height(); i++) {
        for (int j = 0; j < data->Width(); j++) {
            data->QueuePull(i, j);
        }
    }

    data->ProcessPullQueue(bufferPos);

    index = 0;

    for (int i = 0; i < data->Height(); i++) {
        for (int j = 0; j < data->Width(); j++) {
            valuesEigen(i, j) = bufferPos[index];
            index++;
        }
    }

    std::shared_ptr<RBFInterpolation> rbfInterpolator(new RBFInterpolation(rbfFunctionEigen, false, true));
    RBFCoarsening rbfEigen(rbfInterpolator, true, true, false, 1e-3, 0.1, 2, 15, false);

    rbfEigen.compute(positionsEigen, positionsInterpolationEigen);
    rbfEigen.interpolate(valuesEigen, valuesInterpolationEigen);

    std::vector<double> buffer;
    result->ReservePulls(result->Height() * result->Width());

    for (int i = 0; i < result->Height(); i++)
        for (int j = 0; j < result->Width(); j++)
            result->QueuePull(i, j);

    result->ProcessPullQueue(buffer);

    index = 0;

    for (int i = 0; i < result->Height(); i++) {
        for (int j = 0; j < result->Width(); j++) {
            EXPECT_NEAR(valuesInterpolationEigen(i, j), buffer[index], 1e-11);
            index++;
        }
    }
}
