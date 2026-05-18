#include "fks_region.h"

namespace powheg
{
namespace 
{
    bool equalFinalStates(const std::vector<PDG_ID>& fs1, std::vector<PDG_ID> fs2)
    {
        if (fs1.size() != fs2.size())
            return false;

        for (PDG_ID particle : fs1)
        {
            auto it = std::find(fs2.begin(), fs2.end(), particle);
            if (it == fs2.end())
                return false;

            fs2.erase(it);
        }

        return fs2.empty();
    }

    std::optional<ISRCollinearRegion> findLeg1CollinearRegion(
        const BornChannel& bornChannel, 
        RealChannel realChannel, 
        size_t idxEmitted
    )
    {
        PDG_ID emitted = realChannel.outIDs[idxEmitted];

        if (bornChannel.id2 != realChannel.id2)
            return std::nullopt;

        realChannel.outIDs.erase(realChannel.outIDs.begin() + idxEmitted);
        if (!equalFinalStates(bornChannel.outIDs, realChannel.outIDs))
            return std::nullopt;
            
        if (isGluon(emitted) && realChannel.id1 == bornChannel.id1)
        {
            ISRSplitting splitting = isGluon(bornChannel.id1) ? ISRSplitting::GG : ISRSplitting::QQ;
            return ISRCollinearRegion{ splitting, Leg::LEG1 };
        }
        else if (isQuark(emitted))
        {
            if (isGluon(bornChannel.id1) && emitted == realChannel.id1)
                return ISRCollinearRegion{ ISRSplitting::QG, Leg::LEG1 };
            else if (isGluon(realChannel.id1) && isAntiparticle(bornChannel.id1, emitted))
                return ISRCollinearRegion{ ISRSplitting::GQ, Leg::LEG1 };
        }

        return std::nullopt;
    }

    std::optional<ISRCollinearRegion> findLeg2CollinearRegion(
        const BornChannel& bornChannel, 
        RealChannel realChannel, 
        size_t idxEmitted
    )
    {
        PDG_ID emitted = realChannel.outIDs[idxEmitted];

        if (bornChannel.id1 != realChannel.id1)
            return std::nullopt;

        realChannel.outIDs.erase(realChannel.outIDs.begin() + idxEmitted);
        if (!equalFinalStates(bornChannel.outIDs, realChannel.outIDs))
            return std::nullopt;
            
        if (isGluon(emitted) && realChannel.id2 == bornChannel.id2)
        {
            ISRSplitting splitting = isGluon(bornChannel.id2) ? ISRSplitting::GG : ISRSplitting::QQ;
            return ISRCollinearRegion{ splitting, Leg::LEG2 };
        }
        else if (isQuark(emitted))
        {
            if (isGluon(bornChannel.id2) && emitted == realChannel.id2)
                return ISRCollinearRegion{ ISRSplitting::QG, Leg::LEG2 };
            else if (isGluon(realChannel.id2) && isAntiparticle(bornChannel.id2, emitted))
                return ISRCollinearRegion{ ISRSplitting::GQ, Leg::LEG2 };
        }

        return std::nullopt;
    }

} // anonymous namespace

    std::vector<FKSRegion> findRegions(
        const BornChannel& bornChannel,
        const RealChannel& realChannel
    )
    {
        std::vector<FKSRegion> regions;

        for (size_t idxEmitted = 0; idxEmitted < realChannel.outIDs.size(); ++idxEmitted)
        {
            PDG_ID emitted = realChannel.outIDs[idxEmitted];

            if (!isColored(emitted))
                continue;

            auto leg1Coll = findLeg1CollinearRegion(bornChannel, realChannel, idxEmitted);
            if (leg1Coll.has_value())
            {
                regions.push_back({
                    .soft = isGluon(emitted),
                    .isrRegion = leg1Coll,
                    .fsrRegion = std::nullopt
                });
            }

            auto leg2Coll = findLeg2CollinearRegion(bornChannel, realChannel, idxEmitted);
            if (leg2Coll.has_value())
            {
                regions.push_back({
                    .soft = isGluon(emitted),
                    .isrRegion = leg2Coll,
                    .fsrRegion = std::nullopt
                });
            }

            // TODO: Add FSR
        }

        return regions;
    }

} // namespace powheg
