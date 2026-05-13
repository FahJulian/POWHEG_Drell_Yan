#pragma once

#include "powheg_dy/base.h"

namespace powheg_dy
{
    enum class SingularRegionType
    {
        ISR_LEG1,   // emitted parton collinear to incoming leg 1
        ISR_LEG2,   // emitted parton collinear to incoming leg 2
        FSR,        // two collinear outgoing legs
        SOFT,       // emitted gluon soft
    };

    struct SingularRegion
    {
        SingularRegionType type;

        // for ISR: index of the emit
        // for FSR: index of the first collinear leg
        // for SOFT: -1
        int index1;    
        
        // for ISR: index of the radiated parton
        // for FSR: index of the second collinear leg
        // for SOFT: index of the soft parton
        int index2;
    };

    std::vector<SingularRegion> findSingularRegions(
        const BornChannel& born,
        const RealChannel& real
    );

} // namespace powheg_dy
