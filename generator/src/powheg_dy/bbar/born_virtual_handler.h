#pragma once

#include "powheg_dy/base.h"
#include "powheg_dy/config/config.h"
#include "powheg_dy/bbar/bbar_integration_point.h"

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
            const BBarIntegrationPoint& point
        ) const;
        
        double dSigmaVirtual(
            const BBarIntegrationPoint& point, 
            const double muR2
        ) const;

    private:
        const Process& m_process;
        const Config& m_config;
    };
    
} // namespace powheg_dy
