#pragma once

#include "powheg_dy/base.h"
#include "powheg_dy/config/config.h"
#include "powheg_dy/bbar/singular_region.h"
#include "powheg_dy/bbar/bbar_integration_point.h"

namespace powheg_dy
{
    class Process;
    
    class RealMinusCTHandler
    {
    public:
        RealMinusCTHandler(const Process& process, 
            const Config& config)
            : m_process(process), 
                m_config(config)
        {
        }

        double dSigmaRealMinusCT(
            const BBarIntegrationPoint& point,
            const SingularRegion& region
        ) const;

    private:
        const Process& m_process;
        const Config& m_config;
    };
    
} // namespace powheg_dy
