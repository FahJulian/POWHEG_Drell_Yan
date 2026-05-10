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
        void sample(BornPhSpPt& born) const ;

    private:
        double bTilde(const BornPhSpPt& born, const RadiationVariables& rad) const;

    private:
        const Process& m_process;
        const Config& m_config;
        std::shared_ptr<BornPhaseSpace> m_bornPhaseSpace;
        std::shared_ptr<FKSRealPhaseSpace> m_realPhaseSpace;
    };

} // namespace powheg_dy
