#pragma once
 
#include "powheg_dy/phase_space.h"

#include <tuple>
#include <vector>

namespace powheg_dy
{
    class Process;

    struct BornEvent
    {
        int partonId;
        double dSigma;
    };

    class BornEventGenerator
    {
    public:
        BornEventGenerator() = delete;
        BornEventGenerator(const Process& process)
            : m_process(process)
        {
        }

        void computeWeightAndSampleChannel(BornPhaseSpacePt& point) const;

    private:
        std::vector<std::tuple<BornChannel, double>> _computePartonChannelContributions(const BornPhaseSpacePt& point) const;

    private:
        const Process& m_process;
    };

} // namespace powheg_dy
