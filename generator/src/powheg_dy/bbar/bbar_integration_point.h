#pragma once

#include "powheg_dy/base.h"
#include "powheg_dy/phase_space/born_phase_space_point.h"

namespace powheg_dy
{
    // TODO: Split into CollRemnIntPoint and RealIntPoint
    struct BBarIntegrationPoint
    {
        const BornPhSpPt& born;
        const BornChannel& bornChannel;

        double u1 = 0.0;
        double u2 = 0.0;
        double u3 = 0.0;

        double xi = 0.0;
        double y = 0.0;
        double phi = 0.0;

        double unitCubeJacobian = 0.0;

        double zLeg1 = 0.0;
        double zLeg2 = 0.0;

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