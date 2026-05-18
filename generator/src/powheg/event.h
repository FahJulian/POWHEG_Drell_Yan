#pragma once

#include "powheg/base.h"
#include "powheg/emission/emission.h"
#include "powheg/flavour/born_channel.h"
#include "powheg/phase_space/born_phase_space_point.h"
#include "powheg/phase_space/real_phase_space.h"

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

} // namespace powheg
