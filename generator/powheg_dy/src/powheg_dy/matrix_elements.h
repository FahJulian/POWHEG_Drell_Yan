#pragma once

#include "powheg_dy/base.h"

namespace powheg_dy
{
    // TODO: Change to something more general
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

    namespace MatrixElements
    {   
        double virtualFiniteOverBorn(double mSq, double muR2);

        double bbarOverBornApproximation(int partonId, double mSq, double muR2);

    } // namespace MatrixElements

} // namespace powheg_dy
