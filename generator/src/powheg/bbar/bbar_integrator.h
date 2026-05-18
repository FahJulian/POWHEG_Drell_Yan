#pragma once

#include "powheg/config/config.h"
#include "powheg/bbar/born_virtual_handler.h"
#include "powheg/bbar/real_minus_ct_handler.h"
#include "powheg/bbar/bbar_integration_point.h"
#include "powheg/phase_space/real_phase_space.h"
#include "powheg/bbar/collinear_remnants_handler.h"
#include "powheg/phase_space/born_phase_space_point.h"

namespace powheg
{
    class BaseProcess;

    struct BornEvent
    {
        BornPhSpPt phaseSpace;
        BornChannel channel;
        double absoluteWeight;
        int weightSign;
    };

    struct WeightedChannel
    {
        BornChannel channel;
        double amp2Born;
        int weightSign;
        double absoluteWeight;
    };

    class BBarIntegrator
    {
    public: 
        BBarIntegrator(const BaseProcess& process, 
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
            const std::array<double, 3>& unitCube,
            std::vector<WeightedChannel>& weights
        ) const;

    private:
        const BaseProcess& m_process;
        const Config& m_config;
        std::shared_ptr<ISRRealPhaseSpace> m_realPhaseSpace;
        
        const RealMinusCTHandler m_real;
        const BornVirtualHandler m_bornVirtual;
        const CollinearRemnantsHandler m_collRemn;

        int m_sumSigns;
        int m_nEventTrials;
        double m_maxWeight;
    };

} // namespace powheg
