#pragma once

#include "powheg_dy/base.h"
#include "powheg_dy/phase_space/born_phase_space.h"
#include "powheg_dy/phase_space/real_phase_space.h"

namespace powheg_dy
{
    namespace MatrixElements
    {
        double virtualFiniteOverBorn(double mSq, double muR2);

        double bbarOverBornApproximation(int partonId, double mSq, double muR2);

        struct RealOverBornContributions
        {
            double qqbar = 0.0;  // q qbar -> V g
            double gqbar = 0.0;  // g qbar -> V qbar
            double qg    = 0.0;  // q g    -> V q

            double total() const
            {
                return qqbar + gqbar + qg;
            }
        };

        RealOverBornContributions realOverBornContributions(
            const Config& process,
            const RealPhSpPt& real,
            double muF2,
            double muR2,
            bool useCMWALphaS
        );

        double born(const Config& process, const BornPhSpPt& born);

        double realOverBornQQbar(const Config& process, const RealPhSpPt& real, double muR2, bool useCMWALphaS);
        double realOverBornGQbar(const Config& process, const RealPhSpPt& real, double muR2, bool useCMWALphaS);
        double realOverBornQG(const Config& process, const RealPhSpPt& real, double muR2, bool useCMWALphaS);

        double realOverBorn(const Config& process, const RealPhSpPt& real, double muF2, double muR2);
        double realSudakovDensity(const Config& process, const RealPhSpPt& real, double muF2, double muR2);

    } // namespace MatrixElements

} // namespace powheg_dy
