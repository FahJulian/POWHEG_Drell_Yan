#pragma once

#include "powheg/base.h"
#include "powheg/config/config.h"
#include "powheg/bbar/bbar_integration_point.h"

namespace powheg
{
    class BaseProcess;

    class BornVirtualHandler
    {
    public:
        BornVirtualHandler(const BaseProcess& process, 
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
            const double muF2,
            const double muR2
        ) const;

    private:
        const BaseProcess& m_process;
        const Config& m_config;
    };
    
} // namespace powheg
