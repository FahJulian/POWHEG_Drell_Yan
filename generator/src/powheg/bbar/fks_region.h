#pragma once

#include "powheg/base.h"
#include "powheg/flavour/born_channel.h"
#include "powheg/flavour/real_channel.h"
#include "powheg/phase_space/real_phase_space.h"
#include "powheg/phase_space/born_phase_space_point.h"

namespace powheg
{
    struct ISRSingularRegion
    {
        int leg;
        bool isGluon;
    };

    struct FSRSingularRegion
    {
        int emitter;
        int fksParton;
    };

    struct FKSRegionOld
    {
        int fksPartonIdx;
        int emitter;

        bool hasSoft;
        bool hasCollinearLeg1;
        bool hasCollinearLeg2;
        bool hasCollinearFinal;

        // TODO: Directly add flags for ISR splittings here
    };

    std::vector<ISRSingularRegion> findISRSingularRegions(
        const BornChannel& born,
        const RealChannel& real
    );

    std::vector<FSRSingularRegion> findFSRSingularRegions(
        const BornChannel& born,
        const RealChannel& real
    );

} // namespace powheg
