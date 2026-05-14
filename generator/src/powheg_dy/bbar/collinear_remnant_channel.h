#pragma once

#include "powheg_dy/base.h"
#include "powheg_dy/bbar/singular_region.h"
#include "powheg_dy/phase_space/born_phase_space_point.h"
#include "powheg_dy/phase_space/real_phase_space.h"

namespace powheg_dy
{
    enum class CollinearRemnantSplitting
    {
        // q REAL                       qbar BORN                      
        // |                            |
        // |_____g RADIATED     OR      |_____g RADIATED
        // |                            |
        // |                            |
        // q BORN                       qbar REAL
        QQ,       
        
        // q REAL                       g BORN                      
        // |                            |
        // |_____q RADIATED     OR      |_____qbar RADIATED
        // |                            |
        // |                            |
        // g BORN                       qbar REAL
        QG,  
        
        // g REAL                       qbar BORN                      
        // |                            |
        // |_____qbar RADIATED  OR      |_____q RADIATED
        // |                            |
        // |                            |
        // q BORN                       g REAL
        GQ,
        
        // g REAL                       g BORN                      
        // |                            |
        // |_____g    RADIATED  OR      |_____g RADIATED
        // |                            |
        // |                            |
        // g BORN                       g REAL
        GG,
    };

    struct CollinearRemnantChannel
    {
        CollinearRemnantSplitting splitting;
        int leg = 0;
    };

    std::optional<CollinearRemnantChannel> remnantChannelFromRegion(
        const BornChannel& born,
        const RealChannel& real,
        const SingularRegion& region
    );

} // namespace powheg_dy
