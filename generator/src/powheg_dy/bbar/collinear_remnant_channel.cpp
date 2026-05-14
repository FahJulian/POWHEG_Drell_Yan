#include "collinear_remnant_channel.h"

namespace powheg_dy
{
    std::optional<CollinearRemnantChannel> remnantChannelFromRegion(
        const BornChannel& born,
        const RealChannel& real,
        const SingularRegion& region
    )
    {
        // TODO: Implement other regions

        if (region.type == SingularRegionType::ISR_LEG1)
        {
            const int emittedID = real.outIDs[region.index2];

            // q qbar -> X g, i.e. splitting q -> q g on leg 1
            if (real.id1 == born.id1 && emittedID == 21)
                return CollinearRemnantChannel { CollinearRemnantSplitting::QQ, 1 };

            // g qbar -> X qbar, i.e. splitting g -> q qbar on leg 1
            else if (real.id1 == 21 && real.id2 == born.id2 && emittedID == born.id2)
                return CollinearRemnantChannel { CollinearRemnantSplitting::GQ, 1 };

            else 
                throw std::runtime_error("Invalid Singular Region.");
        }
        else if (region.type == SingularRegionType::ISR_LEG2)
        {
            const int emittedID = real.outIDs[region.index2];

            // q qbar -> X g, i.e. splitting q -> q g on leg 2
            if (real.id2 == born.id2 && emittedID == 21)
                return CollinearRemnantChannel { CollinearRemnantSplitting::QQ, 2 };

            // q g -> X q, i.e. splitting g -> qbar q on leg 2
            else if (real.id1 == born.id1 && real.id2 == 21 && emittedID == born.id1)
                return CollinearRemnantChannel { CollinearRemnantSplitting::GQ, 2 };

            else 
                throw std::runtime_error("Invalid Singular Region.");
        }

        return std::nullopt;
    }

} // namespace powheg_dy
