#pragma once

#include "powheg_dy/base.h"
#include "powheg_dy/config.h"
#include "powheg_dy/phase_space/born_phase_space.h"
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
            const BornPhSpPt& born,
            const CollinearRemnantChannel& channel,
            const double u
        ) const;
        
    private:
        double oneMinusZTimesKerelQQ(const double z);
        double oneMinusZTimesKerelQG(const double z);
        double oneMinusZTimesKerelGQ(const double z);
        double oneMinusZTimesKerelGG(const double z);
    
    private:
        const Process& m_process;
        const Config& m_config;
    };
    
} // namespace powheg_dy
