#pragma once

#include "powheg/base.h"
#include "powheg/flavour/pdg_id.h"

namespace powheg
{
    struct RealChannel
    {
        PDG_ID id1;
        PDG_ID id2;
        std::vector<PDG_ID> outIDs;
    };

} // namespace powheg
