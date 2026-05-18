#pragma once

#include "powheg/base.h"
#include "powheg/config/config.h"
#include "powheg/bbar/bbar_cache.h"
#include "powheg/flavour/fks_region.h"
#include "powheg/bbar/collinear_remnant_channel.h"

namespace powheg
{
    class BaseProcess;

    class CollinearRemnantsHandler
    {
    public:
        CollinearRemnantsHandler(const BaseProcess& process, 
            const Config& config)
            : m_process(process), 
                m_config(config)
        {
        }

        double dSigmaCollinearRemnants(
            const BBarCache& cache,
            size_t channelIdx
        ) const;
        
    private:
        double dSigmaCollinearRemnantsLeg1(
            const ISRSplitting& splitting,
            const BBarCache& cache, 
            size_t channelIdx,
            size_t realChannelIdx
        ) const;

        double dSigmaCollinearRemnantsLeg2(
            const ISRSplitting& splitting,
            const BBarCache& cache, 
            size_t channelIdx,
            size_t realChannelIdx
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
        const BaseProcess& m_process;
        const Config& m_config;
    };
    
} // namespace powheg
