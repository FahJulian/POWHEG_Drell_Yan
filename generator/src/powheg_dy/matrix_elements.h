#pragma once

#include "powheg_dy/base.h"

namespace powheg_dy
{
    // TODO: Change to something more general
    struct RealOverBornContributionsOld
    {
        double qqbar = 0.0;
        // Same incoming flavours as the underlying Born event:
        //   q qbar -> V g      or      qbar q -> V g
        
        double gluonLeg1 = 0.0;
        // Gluon replaces incoming Born leg 1:
        //   q qbar -> V     maps to     g qbar -> V qbar
        //   qbar q -> V     maps to     g q    -> V q
        
        double gluonLeg2 = 0.0;
        // Gluon replaces incoming Born leg 2:
        //   q qbar -> V     maps to     q g    -> V q
        //   qbar q -> V     maps to     qbar g -> V qbar

        double total() const
        {
            return qqbar + gluonLeg1 + gluonLeg2;
        }
    };

    namespace MatrixElements
    {   
        double virtualFiniteOverBorn(double mSq, double muR2);

        double bbarOverBornApproximation(int partonId, double mSq, double muR2);

    } // namespace MatrixElements

} // namespace powheg_dy
