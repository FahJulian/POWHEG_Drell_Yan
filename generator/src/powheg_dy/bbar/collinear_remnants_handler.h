#pragma once

#include "powheg_dy/base.h"
#include "powheg_dy/config/config.h"
#include "powheg_dy/bbar/bbar_integration_point.h"
#include "powheg_dy/bbar/collinear_remnant_channel.h"

namespace powheg_dy
{
    class Process;

    class CollinearRemnantsHandler
    {
    public:
        CollinearRemnantsHandler(const Process& process, 
            const Config& config)
            : m_process(process), 
                m_config(config)
        {
        }

        double dSigmaCollinearRemnants(
            const BBarIntegrationPoint& point,
            const CollinearRemnantChannel& channel,
            const double muF2
        ) const;
        
    private:
        double dSigmaCollinearRemnantsLeg1(
            const BBarIntegrationPoint& point,
            const CollinearRemnantSplitting& splitting,
            const double muF2
        ) const;

        double dSigmaCollinearRemnantsLeg2(
            const BBarIntegrationPoint& point,
            const CollinearRemnantSplitting& splitting,
            const double muF2
        ) const;
        
        double oneMinusZTimesPqq(const double z) const;
        double oneMinusZTimesPqg(const double z) const;
        double oneMinusZTimesPgq(const double z) const;
        double oneMinusZTimesPgg(const double z) const;

        double derivativePqq(const double z) const;
        double derivativePqg(const double z) const;
        double derivativePgq(const double z) const;
        double derivativePgg(const double z) const;
    
    private:
        const Process& m_process;
        const Config& m_config;
    };
    
} // namespace powheg_dy
