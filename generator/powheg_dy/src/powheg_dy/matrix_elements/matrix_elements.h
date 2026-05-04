#pragma once

#include "powheg_dy/phase_space/born_phase_space.h"
#include "powheg_dy/phase_space/real_phase_space.h"

namespace powheg_dy
{
    class Process; 

    namespace MatrixElements
    {
        // Angular Born factor without PDFs, flux, phase-space Jacobian, or global \alpha^2 prefactor.
        double bornAngularFactorOld(const Process& process, int partonId, double mSq, double cosTheta);

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

        RealOverBornContributions realOverBornContributionsOld(
            const Process& process,
            const RealPhSpPt& real,
            double muF2,
            double muR2
        );

        double realOverBornQQbarOld(const Process& process, const RealPhSpPt& real, double muR2);
        double realOverBornGQbarOld(const Process& process, const RealPhSpPt& real, double muR2);
        double realOverBornQGOld(const Process& process, const RealPhSpPt& real, double muR2);

        double realOverBornOld(const Process& process, const RealPhSpPt& real, double muF2, double muR2);
        double realSudakovDensityOld(const Process& process, const RealPhSpPt& real, double muF2, double muR2);

        RealOverBornContributions realOverBornContributions(
            const Process& process,
            const RealPhSpPt& real,
            double muF2,
            double muR2
        );

        double born(const Process& process, const BornPhSpPt& real, double muR2);

        double realOverBornQQbar(const Process& process, const RealPhSpPt& real, double muR2);
        double realOverBornGQbar(const Process& process, const RealPhSpPt& real, double muR2);
        double realOverBornQG(const Process& process, const RealPhSpPt& real, double muR2);

        double realOverBorn(const Process& process, const RealPhSpPt& real, double muF2, double muR2);
        double realSudakovDensity(const Process& process, const RealPhSpPt& real, double muF2, double muR2);

    } // namespace MatrixElements

} // namespace powheg_dy
