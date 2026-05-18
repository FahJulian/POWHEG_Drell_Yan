#pragma once

#include "powheg/base.h"
#include "powheg/config/config.h"
#include "powheg/emission/emission.h"
#include "powheg/emission/emission_sampler.h"
#include "powheg/phase_space/born_phase_space_point.h"
#include "powheg/phase_space/real_phase_space.h"

namespace powheg
{
    class BaseProcess;

    struct RealOverBornContributions
    {
        std::vector<std::pair<RealChannel, double>> channels = { };
        double total = 0.0;
    };

    class EmissionGenerator
    {
    public:
        EmissionGenerator(const BaseProcess& process, const Config& config, std::shared_ptr<ISRRealPhaseSpace> realPhaseSpace)
            : m_process(process),
              m_config(config),
              m_sampler(config, realPhaseSpace),
              m_realPhaseSpace(std::move(realPhaseSpace))
        {
        }

        Emission generateEmission(
            const BornPhSpPt& born,
            const BornChannel& bornChannel, 
            const double amp2Born,
            int region
        ) const;
        
    private:
        Emission generateISREmission(
            const BornPhSpPt& born,
            const BornChannel& bornChannel,
            const double amp2Born
        ) const;

        RealOverBornContributions getRealOverBornContributions(
            const RealPhSpPt& real,
            const BornPhSpPt& born,
            const BornChannel& bornChannel,
            const double amp2Born,
            const double muF2,
            const double muR2
        ) const;
        
    private:
        const BaseProcess& m_process;
        const Config& m_config;
        const EmissionSampler m_sampler;
        std::shared_ptr<ISRRealPhaseSpace> m_realPhaseSpace;
    };

} // namespace powheg
