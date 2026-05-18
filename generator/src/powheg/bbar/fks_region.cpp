#include "fks_region.h"

#include "powheg/phase_space/born_phase_space_point.h"
#include "powheg/phase_space/real_phase_space.h"

namespace powheg
{
    std::vector<ISRSingularRegion> findISRSingularRegions(
        const BornChannel& born,
        const RealChannel& real
    )
    {
        std::vector<ISRSingularRegion> regions;

        // for (int outID : real.outIDs)
        // {
        //     if (outID == 21)
        //     {
        //         if (real.id1)
        //     }
        // }
    }

    std::vector<FSRSingularRegion> findFSRSingularRegions(
        const BornChannel& born,
        const RealChannel& real
    )
    {
        // TODO: Implement
        return { };
    }

} // namespace powheg
