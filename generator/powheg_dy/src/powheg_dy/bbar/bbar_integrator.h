#pragma once

#include "powheg_dy/config.h"
#include "powheg_dy/phase_space/born_phase_space.h"
#include "powheg_dy/phase_space/real_phase_space.h"

namespace powheg_dy
{
    class Process;

    class BBarIntegrator
    {
    public: 
        BBarIntegrator(const Process& process, 
            const Config& config, 
            std::shared_ptr<BornPhaseSpace> bornPhaseSpace, 
            std::shared_ptr<FKSRealPhaseSpace> realPhaseSpace)
            : m_process(process), 
                m_config(config), 
                m_bornPhaseSpace(std::move(bornPhaseSpace)), 
                m_realPhaseSpace(std::move(realPhaseSpace))
        {
        }

        void clear();

        void determineMaxWeight();
        void computeWeightAndSampleChannel(BornPhSpPt& born) const ;
        BornPhSpPt sampleAccordingtoBTilde();

        double getAcceptanceRatio() const;
        double getTotalCrossSection() const;

    private:
        double bTilde(const BornPhSpPt& born, const RadiationVariables& rad) const;

    private:
        const Process& m_process;
        const Config& m_config;
        std::shared_ptr<BornPhaseSpace> m_bornPhaseSpace;
        std::shared_ptr<FKSRealPhaseSpace> m_realPhaseSpace;

        int m_nEventTrials;
        double m_maxWeight;
    };

} // namespace powheg_dy
