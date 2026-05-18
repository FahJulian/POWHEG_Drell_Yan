#pragma once

#include "powheg/base.h"
#include "powheg/config/config.h"
#include "powheg/bbar/bbar_cache.h"

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
            const BBarCache& cache,
            size_t channelIdx
        ) const;
        
        double dSigmaVirtual(
            const BBarCache& cache,
            size_t channelIdx
        ) const;

    private:
        const BaseProcess& m_process;
        const Config& m_config;
    };
    
} // namespace powheg
