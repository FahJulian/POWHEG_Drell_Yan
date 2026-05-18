#include "channel.h"

namespace powheg
{
namespace 
{
    bool hasCollinearLeg1(const std::vector<FKSRegion>& regions)
    {
        for (const FKSRegion& region : regions)
        {
            if (region.isrRegion.has_value() && region.isrRegion->leg == Leg::LEG1)
                return true;
        }

        return false;
    }

    bool hasCollinearLeg2(const std::vector<FKSRegion>& regions)
    {
        for (const FKSRegion& region : regions)
        {
            if (region.isrRegion.has_value() && region.isrRegion->leg == Leg::LEG2)
                return true;
        }

        return false;
    }

} // anonymous namespace

    std::vector<Channel> findChannels(
        const std::vector<BornChannel>& bornChannels,
        const std::vector<RealChannel>& realChannels
    )
    {
        std::vector<Channel> channels;

        for (const BornChannel& bornChannel : bornChannels)
        {
            Channel channel{ bornChannel };

            for (const RealChannel& realChannel : realChannels)
            {
                auto regions = findRegions(bornChannel, realChannel);
                if (!regions.empty())
                {
                    channel.fksRegions.push_back(regions);
                    channel.relevantRealChannels.push_back(realChannel);
                    channel.hasCollinearLeg1.push_back(hasCollinearLeg1(regions));
                    channel.hasCollinearLeg2.push_back(hasCollinearLeg2(regions));
                }
            }

            channels.push_back(channel);
        }

        return channels;
    }

} // namespace powheg
