#pragma once

#include "powheg/base.h"
#include "powheg/flavour/pdg_id.h"
#include "powheg/flavour/born_channel.h"
#include "powheg/flavour/real_channel.h"

namespace powheg
{
    enum class ISRSplitting
    {
        QQ,
        QG,
        GQ,
        GG
    };

    enum class Leg 
    {
        LEG1,
        LEG2
    };

    struct ISRCollinearRegion
    {
        ISRSplitting splitting;
        Leg leg;
    };

    struct FSRCollinearRegion
    {
        size_t idxEmitter = -1;
        size_t idxEmitted = -1;
    };

    struct FKSRegion
    {
        bool soft = false;
        std::optional<ISRCollinearRegion> isrRegion = std::nullopt;
        std::optional<FSRCollinearRegion> fsrRegion = std::nullopt;
    };

    std::vector<FKSRegion> findRegions(const BornChannel& bornChannel, const RealChannel& realChannel);

} // namespace powheg
