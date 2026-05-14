#pragma once

#include "powheg_dy/config/config.h"
#include "powheg_dy/bbar/born_virtual_handler.h"
#include "powheg_dy/bbar/real_minus_ct_handler.h"
#include "powheg_dy/bbar/bbar_integration_point.h"
#include "powheg_dy/phase_space/real_phase_space.h"
#include "powheg_dy/bbar/collinear_remnants_handler.h"
#include "powheg_dy/phase_space/born_phase_space_point.h"

namespace powheg_dy
{
    class Process;

    struct BornEvent
    {
        BornPhSpPt phaseSpace;
        BornChannel channel;
        double amp2Born;
        double absoluteWeight;
        int weightSign;
    };

    class BBarIntegrator
    {
    public: 
        BBarIntegrator(const Process& process, 
            const Config& config, 
            std::shared_ptr<ISRRealPhaseSpace> realPhaseSpace)
            : m_process(process), 
                m_config(config), 
                m_realPhaseSpace(std::move(realPhaseSpace)),
                m_real(process, config),
                m_bornVirtual(process, config),
                m_collRemn(process, config)
        {
        }

        void clear();

        void determineMaxWeight();
        BornEvent sampleAccordingtoBTilde();

        double getAcceptanceRatio() const;
        double getTotalCrossSection() const;
        double getAbsCrossSection() const;

    private:
        BornEvent computeWeightAndSampleChannel(
            const BornPhSpPt& born
        ) const;

        BBarIntegrationPoint generateIntegrationPoint(
            const BornPhSpPt& born, 
            const BornChannel& bornChannel, 
            const double amp2Born,
            const double muF2,
            const double muR2,
            const std::array<double, 3>& unitCube
        ) const;

        double bTilde(
            const BornPhSpPt& born,
            const double amp2Born, 
            const BornChannel& bornChannel, 
            const std::array<double, 3>& unitCube
        ) const;

    private:
        const Process& m_process;
        const Config& m_config;
        std::shared_ptr<ISRRealPhaseSpace> m_realPhaseSpace;
        
        const RealMinusCTHandler m_real;
        const BornVirtualHandler m_bornVirtual;
        const CollinearRemnantsHandler m_collRemn;

        int m_sumSigns;
        int m_nEventTrials;
        double m_maxWeight;
    };

} // namespace powheg_dy
