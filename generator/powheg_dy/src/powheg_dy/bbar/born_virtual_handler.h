#pragma once

#include "powheg_dy/base.h"
#include "powheg_dy/config.h"
#include "powheg_dy/phase_space/born_phase_space.h"

namespace powheg_dy
{
    class Process;

    class BornVirtualHandler
    {
    public:
        BornVirtualHandler(const Process& process, 
            const Config& config)
            : m_process(process), 
                m_config(config)
        {
        }

        double dSigmaBorn(
            const BornPhSpPt& born
        ) const;
        
        double dSigmaVirtual(
            const BornPhSpPt& born, 
            const double muR2
        ) const;

    private:
        const Process& m_process;
        const Config& m_config;
    };
    
} // namespace powheg_dy
