#pragma once

#include "powheg/base.h"
#include "powheg/flavour/fks_region.h"
#include "powheg/flavour/born_channel.h"
#include "powheg/flavour/real_channel.h"

namespace powheg
{
    struct Channel
    {
        BornChannel bornChannel;
        std::vector<RealChannel> relevantRealChannels;
        std::vector<std::vector<FKSRegion>> fksRegions;

        std::vector<bool> hasCollinearLeg1;
        std::vector<bool> hasCollinearLeg2;
    };

    std::vector<Channel> findChannels(
        const std::vector<BornChannel>& bornChannels,
        const std::vector<RealChannel>& realChannels
    );

} // namespace powheg
