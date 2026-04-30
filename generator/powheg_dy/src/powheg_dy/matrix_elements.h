#pragma once

#include "phase_space.h"

namespace powheg_dy
{
    class Process; 

    namespace MatrixElements
    {
        // Angular Born factor without PDFs, flux, phase-space Jacobian, or global alpha^2 prefactor.
        double bornAngularFactor(const Process& process, int partonId, double mSq, double cosTheta);

        // Finite Drell-Yan virtual correction divided by the Born term, in the same convention
        // as POWHEG's virtual.f: the global alpha_s/(2*pi) factor is NOT included here.
        double virtualFiniteOverBorn(double mSq, double muR2);

        // Placeholder for the future true POWHEG Bbar/B. For now this deliberately returns 1.
        // Do not mistake this for a full NLO K-factor: the integrated real/subtraction and
        // collinear-remnant terms are not implemented yet.
        double bbarOverBornApproximation(int partonId, double mSq, double muR2);

    } // namespace MatrixElements

} // namespace powheg_dy
