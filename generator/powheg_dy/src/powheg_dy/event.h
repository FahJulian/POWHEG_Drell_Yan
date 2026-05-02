#pragma once

#include "powheg_dy/emission.h"
#include "powheg_dy/phase_space.h"
#include "powheg_dy/real_phase_space.h"

namespace powheg_dy
{
    struct Event
    {
        BornPhSpPt born;
        RealPhSpPt real;
        Emission emission;
    };

} // namespace powheg_dy
