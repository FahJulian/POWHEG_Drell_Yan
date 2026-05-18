#pragma once

#include "powheg_dy/base.h"
#include "powheg_dy/phase_space/real_phase_space.h"
#include "powheg_dy/phase_space/born_phase_space_point.h"

namespace powheg
{
    // TODO: Split into CollRemnIntPoint and RealIntPoint
    struct BBarIntegrationPoint
    {
        BornPhSpPt born;
        BornChannel bornChannel;
        RealPhSpPt real;

        double u1 = 0.0;
        double u2 = 0.0;
        double u3 = 0.0;

        double xi = 0.0;
        double y = 0.0;
        double phi = 0.0;

        double xiTilde = 0.0;
        double xiLeg1 = 0.0;
        double xiLeg2 = 0.0;

        double xiMax = 0.0;
        double xiMaxLeg1 = 0.0;
        double xiMaxLeg2 = 0.0;

        double jacobian = 0.0;
        double jacobianOverXiSqSoft = 0.0;
        double jacobianLeg1 = 0.0;
        double jacobianLeg2 = 0.0;

        double zLeg1CollRemn = 0.0;
        double zLeg2CollRemn = 0.0;

        double f1Born = 0.0;
        double f2Born = 0.0;

        double f1RealQ = 0.0;
        double f2RealQ = 0.0;

        double f1RealG = 0.0;
        double f2RealG = 0.0;

        double alphaS = 0.0;

        double amp2Born = 0.0;
    };

} // namespace powheg_dy