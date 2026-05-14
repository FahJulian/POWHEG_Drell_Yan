#include "singular_region.h"

#include "powheg_dy/phase_space/born_phase_space_point.h"
#include "powheg_dy/phase_space/real_phase_space.h"

namespace powheg_dy
{
    std::vector<SingularRegion> findSingularRegions(
        const BornChannel& born,
        const RealChannel& real
    )
    {
        std::vector<SingularRegion> regions;

        // TODO: Add FSR

        // q qbar -> X g
        if (real.id1 == born.id1
            && real.id2 == born.id2
            && real.outIDs.size() == born.outIDs.size() + 1
            && real.outIDs.back() == 21)
        {
            regions.push_back({
                .type = SingularRegionType::SOFT,
                .index1 = -1,
                .index2 = static_cast<int>(real.outIDs.size()) - 1,
            });

            regions.push_back({
                .type = SingularRegionType::ISR_LEG1,
                .index1 = 0,
                .index2 = static_cast<int>(real.outIDs.size()) - 1,
            });

            regions.push_back({
                .type = SingularRegionType::ISR_LEG2,
                .index1 = 1,
                .index2 = static_cast<int>(real.outIDs.size()) - 1,
            });

            return regions;
        }

        // q g -> X q
        if (real.id1 == born.id1
            && real.id2 == 21
            && real.outIDs.size() == born.outIDs.size() + 1
            && real.outIDs.back() == born.id1)
        {
            regions.push_back({
                .type = SingularRegionType::ISR_LEG2,
                .index1 = 1,
                .index2 = static_cast<int>(real.outIDs.size()) - 1,
            });

            return regions;
        }

        // g qbar -> X qbar
        if (real.id1 == 21
            && real.id2 == born.id2
            && real.outIDs.size() == born.outIDs.size() + 1
            && real.outIDs.back() == born.id2)
        {
            regions.push_back({
                .type = SingularRegionType::ISR_LEG1,
                .index1 = 0,
                .index2 = static_cast<int>(real.outIDs.size()) - 1,
            });

            return regions;
        }

        // TODO: Add FSR

        throw std::runtime_error("Real channel does not match expected DY real channels.");
    }

} // namespace powheg_dy
