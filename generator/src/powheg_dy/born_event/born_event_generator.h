#pragma once
 
#include "powheg_dy/base.h"
#include "powheg_dy/config.h"
#include "powheg_dy/phase_space/born_phase_space.h"

namespace powheg_dy
{
    class Process;

    class BornEventGenerator
    {
    public:
        BornEventGenerator() = delete;
        BornEventGenerator(const Process& process, const Config& config, std::shared_ptr<BornPhaseSpace> phaseSpace)
            : m_process(process), m_config(config), m_bornPhaseSpace(std::move(phaseSpace))
        {
        }

        void computeWeightAndSampleChannel(BornPhSpPt& point) const;

    private:
        std::vector<std::tuple<BornPhSpPt, double>> computePartonChannelContributions(const BornPhSpPt& point) const;

    private:
        const Process& m_process;
        const Config& m_config;
        std::shared_ptr<BornPhaseSpace> m_bornPhaseSpace;
    };

} // namespace powheg_dy
