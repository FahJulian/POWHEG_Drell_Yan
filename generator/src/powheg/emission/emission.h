#pragma once

#include "powheg/base.h"
#include "powheg/flavour/real_channel.h"
#include "powheg/phase_space/real_phase_space.h"

namespace powheg
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

} // namespace powheg
