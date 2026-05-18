#pragma once

#include "powheg/base.h"
#include "powheg/config/config.h"
#include "powheg/bbar/bbar_cache.h"
#include "powheg/flavour/fks_region.h"

namespace powheg
{
    class BaseProcess;
    
    class RealMinusCTHandler
    {
    public:
        RealMinusCTHandler(
            const BaseProcess& process, 
            const Config& config, 
            const ISRRealPhaseSpace& realPhaseSpace)
            : m_process(process), 
                m_config(config),
                m_realPhaseSpace(realPhaseSpace)
        {
        }

        double dSigmaRealMinusCT(
            const BBarCache& cache,
            size_t channelIdx
        ) const;

    private:
        double dSigmaReal(
            const BBarCacheOld& cache,
            const RealChannel& realChannel,
            const double muF2
        ) const;

        double softCounterterm(
            const BBarCache& cache,
            size_t channelIdx,
            const FKSRegion& region
        ) const;

        double leg1CollinearCounterterm(
            const BBarCache& cache,
            size_t channelIdx,
            size_t realChannelIdx,
            const FKSRegion& region
        ) const;

        double leg2CollinearCounterterm(
            const BBarCache& cache,
            size_t channelIdx,
            size_t realChannelIdx,
            const FKSRegion& region
        ) const;

        double leg1SoftCollinearCounterterm(
            const BBarCache& cache,
            size_t channelIdx,
            const FKSRegion& region
        ) const;

        double leg2SoftCollinearCounterterm(
            const BBarCache& cache,
            size_t channelIdx,
            const FKSRegion& region
        ) const;

        double fksPartition(
            const FKSRegion& region,
            double y
        ) const;

        double oneMinusZTimesPqq(const double z) const;
        double oneMinusZTimesPqg(const double z) const;
        double oneMinusZTimesPgq(const double z) const;
        double oneMinusZTimesPgg(const double z) const;

        double oneMinusZTimesPqqAtOne() const;
        double oneMinusZTimesPqgAtOne() const;
        double oneMinusZTimesPgqAtOne() const;
        double oneMinusZTimesPggAtOne() const;

        const BaseProcess& m_process;
        const Config& m_config;
        const ISRRealPhaseSpace& m_realPhaseSpace;
    };
    
} // namespace powheg
