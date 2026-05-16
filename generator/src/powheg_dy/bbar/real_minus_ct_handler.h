#pragma once

#include "powheg_dy/base.h"
#include "powheg_dy/config/config.h"
#include "powheg_dy/bbar/fks_region.h"
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
            const RealChannel& realChannel,
            const FKSRegion& region,
            const double muF2
        ) const;

    private:
        double dSigmaReal(
            const BBarIntegrationPoint& point,
            const RealChannel& realChannel,
            const double muF2
        ) const;

        double softCounterterm(
            const BBarIntegrationPoint& point,
            const RealChannel& realChannel,
            const FKSRegion& region
        ) const;

        double leg1CollinearCounterterm(
            const BBarIntegrationPoint& point,
            const RealChannel& realChannel,
            const FKSRegion& region,
            const double muF2
        ) const;

        double leg2CollinearCounterterm(
            const BBarIntegrationPoint& point,
            const RealChannel& realChannel,
            const FKSRegion& region,
            const double muF2
        ) const;

        double leg1SoftCollinearCounterterm(
            const BBarIntegrationPoint& point,
            const RealChannel& realChannel,
            const FKSRegion& region
        ) const;

        double leg2SoftCollinearCounterterm(
            const BBarIntegrationPoint& point,
            const RealChannel& realChannel,
            const FKSRegion& region
        ) const;

        double oneMinusZTimesPqq(const double z) const;
        double oneMinusZTimesPqg(const double z) const;
        double oneMinusZTimesPgq(const double z) const;
        double oneMinusZTimesPgg(const double z) const;

        double oneMinusZTimesPqqAtOne() const;
        double oneMinusZTimesPqgAtOne() const;
        double oneMinusZTimesPgqAtOne() const;
        double oneMinusZTimesPggAtOne() const;

        const Process& m_process;
        const Config& m_config;
    };
    
} // namespace powheg_dy
