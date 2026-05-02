#pragma once

#include "real_phase_space.h"

namespace powheg_dy
{
    enum class RadiationRegion
    {
        None,

        // Initial-state singular direction.
        // Plus  means y -> +1, i.e. emitted parton collinear to incoming leg 1.
        // Minus means y -> -1, i.e. emitted parton collinear to incoming leg 2.
        Plus,
        Minus
    };

    struct Emission
    {
        RadiationVariables rad;

        RadiationRegion reg = RadiationRegion::None;
        double kt2 = 0.0;

        double exactDensity = 0.0;
        double upperDensity = 0.0;
        double acceptanceRatio = 0.0;

        bool rejected = false;

        Emission reject()
        {
            rad = {};
            reg = RadiationRegion::None;
            kt2 = 0.0;
            exactDensity = 0.0;
            upperDensity = 0.0;
            acceptanceRatio = 0.0;
            rejected = true;
            return *this;
        }
    };

} // namespace powheg_dy
