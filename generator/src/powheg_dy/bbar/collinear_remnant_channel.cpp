#include "collinear_remnant_channel.h"

namespace powheg_dy
{
    std::vector<CollinearRemnantChannel> remnantChannelsFromRegion(
        const BornChannel& born,
        const RealChannel& real,
        const FKSRegion& region
    )
    {
        // TODO: Implement other regions

        std::vector<CollinearRemnantChannel> channels;

        if (region.hasCollinearLeg1)
        {
            const int emittedID = real.outIDs[region.fksPartonIdx];

            // q qbar -> X g, i.e. splitting q -> q g on leg 1
            if (real.id1 == born.id1 && real.id2 == born.id2 && emittedID == 21)
                channels.push_back({ CollinearRemnantSplitting::QQ, 1 });

            // g qbar -> X qbar, i.e. splitting g -> q qbar on leg 1
            else if (real.id1 == 21 && real.id2 == born.id2 && emittedID == -born.id1)
                channels.push_back({ CollinearRemnantSplitting::GQ, 1 });

            else 
                throw std::runtime_error("Invalid FKS Region.");
        }

        if (region.hasCollinearLeg2)
        {
            const int emittedID = real.outIDs[region.fksPartonIdx];

            // q qbar -> X g, i.e. splitting q -> q g on leg 2
            if (real.id1 == born.id1 && real.id2 == born.id2 && emittedID == 21)
                channels.push_back({ CollinearRemnantSplitting::QQ, 2 });

            // q g -> X q, i.e. splitting g -> qbar q on leg 2
            else if (real.id1 == born.id1 && real.id2 == 21 && emittedID == -born.id2)
                channels.push_back({ CollinearRemnantSplitting::GQ, 2 });

            else 
                throw std::runtime_error("Invalid FKS Region.");
        }

        return channels;
    }

} // namespace powheg_dy
