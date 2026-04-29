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

        BornEvent computeWeightAndSampleParton(const PhaseSpacePoint& point) const;

    private:
        std::vector<std::tuple<int, double>> _computePartonChannelContributions(const PhaseSpacePoint& point) const;
        std::tuple<double, double> _neutralCurrentCouplingFactors(bool upType, double mSq) const;

    private:
        const Process& m_process;
    };

} // namespace powheg_dy
