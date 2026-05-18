#pragma once

#include "powheg_dy/base.h"
#include "powheg_dy/emission/emission.h"
#include "powheg_dy/phase_space/born_phase_space_point.h"
#include "powheg_dy/phase_space/real_phase_space.h"

namespace powheg
{
    struct Event
    {
        BornPhSpPt born;
        BornChannel bornChannel;
        int weightSign;
        
        RealPhSpPt real;
        Emission emission;
    };

} // namespace powheg_dy
