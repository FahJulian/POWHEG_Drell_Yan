#pragma once

#include "real_phase_space.h"

namespace powheg_dy
{
    struct RadiationChannel
    {
        int idPlus = 0;
        int idMinus = 0;
        int idRadiated = 0;
    };

    struct Emission
    {
        RadiationVariables rad;
        RadiationChannel channel;

        double kt2 = 0.0;

        double exactDensity = 0.0;
        double upperDensity = 0.0;
        double acceptanceRatio = 0.0;

        bool rejected = false;

        Emission reject()
        {
            rad = {};
            channel = { };
            kt2 = 0.0;
            exactDensity = 0.0;
            upperDensity = 0.0;
            acceptanceRatio = 0.0;
            rejected = true;
            return *this;
        }
    };

} // namespace powheg_dy
