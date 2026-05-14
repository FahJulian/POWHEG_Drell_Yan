#pragma once

#include "powheg_dy/base.h"
#include "powheg_dy/phase_space/real_phase_space.h"

namespace powheg_dy
{
    struct RadiationChannel
    {
        int id1 = 0;
        int id2 = 0;
        int idRadiated = 0;
    };

    struct Emission
    {
        RadiationVariables rad;
        RealChannel channel;

        double kt2 = 0.0;
        bool rejected = false;

        Emission reject()
        {
            rad = {};
            channel = { };
            kt2 = 0.0;
            rejected = true;
            return *this;
        }
    };

} // namespace powheg_dy
