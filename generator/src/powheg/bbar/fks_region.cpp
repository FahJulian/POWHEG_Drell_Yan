#include "fks_region.h"

#include "powheg/phase_space/born_phase_space_point.h"
#include "powheg/phase_space/real_phase_space.h"

namespace powheg
{
    std::vector<FKSRegion> findSingularRegions(
        const BornChannel& born,
        const RealChannel& real
    )
    {
        throw std::runtime_error("Not yet implemented.");
        return { };
    }

} // namespace powheg
