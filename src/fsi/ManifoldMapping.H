
/*
 * Copyright [2016] <David Blom>
 */

#ifndef ManifoldMapping_H
#define ManifoldMapping_H

#include "SpaceMapping.H"
#include "ImplicitMultiLevelFsiSolver.H"

namespace fsi
{
    class ManifoldMapping : public SpaceMapping
    {
        public:
            using SpaceMapping::performPostProcessing;

            ManifoldMapping(
                shared_ptr<SurrogateModel> fineModel,
                shared_ptr<SurrogateModel> coarseModel,
                int maxIter
                );

            ManifoldMapping(
                shared_ptr<SurrogateModel> fineModel,
                shared_ptr<SurrogateModel> coarseModel,
                int maxIter,
                int maxUsedIterations,
                int nbReuse,
                int reuseInformationStartingFromTimeIndex,
                scalar singularityLimit,
                bool updateJacobian
                );

            ManifoldMapping(
                shared_ptr<SurrogateModel> fineModel,
                shared_ptr<SurrogateModel> coarseModel,
                int maxIter,
                int maxUsedIterations,
                int nbReuse,
                int reuseInformationStartingFromTimeIndex,
                scalar singularityLimit,
                bool updateJacobian,
                bool initialSolutionCoarseModel
                );

            virtual ~ManifoldMapping();

            virtual void performPostProcessing(
                const vector & y,
                const vector & x0,
                vector & xk,
                bool residualCriterium
                );

            void applyScaling( vector & vec );

            void applyScaling( matrix & mat );

            void determineScalingFactors( const vector & output );

            const bool updateJacobian;
            const bool initialSolutionCoarseModel;
            bool scaling;
            int iter;
            vector scalingFactors;
            matrix Tkprev;
            int sizeVar0;
            int sizeVar1;
    };
} // namespace fsi

#endif
