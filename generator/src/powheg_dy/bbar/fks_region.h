#pragma once

#include "powheg_dy/base.h"
#include "powheg_dy/phase_space/born_phase_space_point.h"
#include "powheg_dy/phase_space/real_phase_space.h"

namespace powheg_dy
{
    struct FKSRegion
    {
        int fksPartonIdx;
        int emitter;

        bool hasSoft;
        bool hasCollinearLeg1;
        bool hasCollinearLeg2;
        bool hasCollinearFinal;

        // TODO: Directly add flags for ISR splittings here
    };

    std::vector<FKSRegion> findSingularRegions(
        const BornChannel& born,
        const RealChannel& real
    );

} // namespace powheg_dy
