#pragma once

#include "powheg_dy/phase_space.h"
#include "powheg_dy/real_phase_space.h"

namespace powheg_dy
{
    class Process; 

    namespace MatrixElements
    {
        // Angular Born factor without PDFs, flux, phase-space Jacobian, or global \alpha^2 prefactor.
        double bornAngularFactor(const Process& process, int partonId, double mSq, double cosTheta);

        // Finite Drell-Yan virtual correction divided by the Born term, in the same convention
        // as POWHEG's virtual.f: the global \alpha_s/(2\pi) factor is NOT included here.
        double virtualFiniteOverBorn(double mSq, double muR2);

        // Placeholder for the future true POWHEG Bbar/B. For now this deliberately returns 1.
        // Do not mistake this for a full NLO K-factor: the integrated real/subtraction and
        // collinear-remnant terms are not implemented yet.
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
            const Process& process,
            const RealPhSpPt& real,
            double muF2,
            double muR2
        );

        double realOverBornQQbar(const Process& process, const RealPhSpPt& real, double muR2);
        double realOverBornGQbar(const Process& process, const RealPhSpPt& real, double muR2);
        double realOverBornQG(const Process& process, const RealPhSpPt& real, double muR2);

        double realOverBorn(const Process& process, const RealPhSpPt& real, double muF2, double muR2);
        double realSudakovDensity(const Process& process, const RealPhSpPt& real, double muF2, double muR2);

    } // namespace MatrixElements

} // namespace powheg_dy
