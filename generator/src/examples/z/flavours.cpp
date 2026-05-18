#include "drell_yan.h"

namespace powheg
{
    std::vector<BornChannel> DrellYanProcess::bornChannels() const
    {
        return {
            { PDG_ID::D, PDG_ID::DBAR, { PDG_ID::E_MINUS, PDG_ID::E_PLUS } },
            { PDG_ID::U, PDG_ID::UBAR, { PDG_ID::E_MINUS, PDG_ID::E_PLUS } },
            { PDG_ID::S, PDG_ID::SBAR, { PDG_ID::E_MINUS, PDG_ID::E_PLUS } },
            { PDG_ID::C, PDG_ID::CBAR, { PDG_ID::E_MINUS, PDG_ID::E_PLUS } },
            { PDG_ID::B, PDG_ID::BBAR, { PDG_ID::E_MINUS, PDG_ID::E_PLUS } },
            { PDG_ID::DBAR, PDG_ID::D, { PDG_ID::E_MINUS, PDG_ID::E_PLUS } },
            { PDG_ID::UBAR, PDG_ID::U, { PDG_ID::E_MINUS, PDG_ID::E_PLUS } },
            { PDG_ID::SBAR, PDG_ID::S, { PDG_ID::E_MINUS, PDG_ID::E_PLUS } },
            { PDG_ID::CBAR, PDG_ID::C, { PDG_ID::E_MINUS, PDG_ID::E_PLUS } },
            { PDG_ID::BBAR, PDG_ID::B, { PDG_ID::E_MINUS, PDG_ID::E_PLUS } }
        };
    }

    std::vector<RealChannel> DrellYanProcess::realChannels() const
    {   
        std::vector<RealChannel> channels;

        for (const BornChannel& bornChannel : bornChannels())
        {
            const PDG_ID outID1 = bornChannel.outIDs[0];
            const PDG_ID outID2 = bornChannel.outIDs[1];

            channels.push_back({ bornChannel.id1, bornChannel.id2, { outID1, outID2, PDG_ID::G } });
            channels.push_back({ bornChannel.id1, PDG_ID::G,       { outID1, outID2, bornChannel.id1 } });
            channels.push_back({ PDG_ID::G,       bornChannel.id2, { outID1, outID2, bornChannel.id2 } });
        }

        return channels;
    }

    std::vector<RealChannel> DrellYanProcess::realChannelsOld(const BornChannel& bornChannel) const
    {   
        const PDG_ID outID1 = bornChannel.outIDs[0];
        const PDG_ID outID2 = bornChannel.outIDs[1];

        return {
            { bornChannel.id1, bornChannel.id2, { outID1, outID2, PDG_ID::G } },
            { bornChannel.id1, PDG_ID::G,       { outID1, outID2, bornChannel.id1 } },
            { PDG_ID::G,       bornChannel.id2, { outID1, outID2, bornChannel.id2 } }
        };
    }

} // namespace powheg
