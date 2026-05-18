#pragma once

#include "powheg/config/config.h"
#include "powheg/bbar/bbar_cache.h"
#include "powheg/flavour/fks_region.h"
#include "powheg/bbar/born_virtual_handler.h"
#include "powheg/bbar/real_minus_ct_handler.h"
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
            const ISRRealPhaseSpace& realPhaseSpace)
            : m_process(process), 
                m_config(config), 
                m_realPhaseSpace(realPhaseSpace),
                m_real(process, config, realPhaseSpace),
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
        BBarCache prepareCache(
            const std::array<double, 3>& unitCube,
            const BornPhSpPt& born
        ) const;

        ChannelCache prepareChannelCache(
            const BBarCache& cache, 
            const Channel& channel
        ) const;

        void prepareCollinearLuminosities(
            ChannelCache& cache,
            const BornPhSpPt& born, 
            const Channel& channel,
            double muF2
        ) const;

        RadiationCache sampleRadiationVariables(
            const std::array<double, 3>& unitCube,
            const BornPhSpPt& born
        ) const;

        BornEvent computeWeightAndSampleChannel(
            const BornPhSpPt& born
        ) const;

        BBarCacheOld prepareCache(
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
        const ISRRealPhaseSpace& m_realPhaseSpace;
        
        const RealMinusCTHandler m_real;
        const BornVirtualHandler m_bornVirtual;
        const CollinearRemnantsHandler m_collRemn;

        int m_sumSigns;
        int m_nEventTrials;
        double m_maxWeight;
    };

} // namespace powheg
