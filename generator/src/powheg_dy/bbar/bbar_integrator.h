#pragma once

#include "powheg_dy/config/config.h"
#include "powheg_dy/bbar/born_virtual_handler.h"
#include "powheg_dy/bbar/collinear_remnants_handler.h"
#include "powheg_dy/bbar/bbar_integration_point.h"
#include "powheg_dy/phase_space/born_phase_space_point.h"
#include "powheg_dy/phase_space/real_phase_space.h"

namespace powheg_dy
{
    class Process;

    struct BornEvent
    {
        BornPhSpPt phaseSpace;
        BornChannel channel;
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
        BornEvent computeWeightAndSampleChannel(const BornPhSpPt& born) const;

        double bTildeOld(BornPhSpPt& born, const std::array<double, 3>& unitCube, const std::array<double, 4>& unitX) const;

        double bornPlusVirtualContribution(const BornPhSpPt& born, double muF2, double muR2) const;
        double realMinusCounterTermContribution(const BornPhSpPt& born, const std::array<double, 3>& unitCube, double muF2, double muR2) const;

        double luminosity(const BornPhSpPt& born, const RadiationVariables& rad, int id1, int id2, double muF2) const;

        double RHatqqbarPlus(const BornPhSpPt& born, const std::array<double, 3>& unitCube, double muF2, double muR2) const;
        double RHatqqbarMinus(const BornPhSpPt& born, const std::array<double, 3>& unitCube, double muF2, double muR2) const;

        double AqqbarPlusLimitAware(const BornPhSpPt& born, const RadiationVariables& rad, double muR2) const;
        double AqqbarMinusLimitAware(const BornPhSpPt& born, const RadiationVariables& rad, double muR2) const;
        
        double AqqbarPlus(const RealPhSpPt& real, double muR2) const;
        double AqqbarMinus(const RealPhSpPt& real, double muR2) const;

        double AqqbarPlusCollinearLimit(const BornPhSpPt& born, const RadiationVariables& rad, double muR2) const;
        double AqqbarMinusCollinearLimit(const BornPhSpPt& born, const RadiationVariables& rad, double muR2) const;

        double AqqbarPlusSoftLimit(const BornPhSpPt& born, const RadiationVariables& rad, double muR2) const;
        double AqqbarMinusSoftLimit(const BornPhSpPt& born, const RadiationVariables& rad, double muR2) const;

        double AqqbarPlusSoftCollinearLimit(const BornPhSpPt& born, double muR2) const;
        double AqqbarMinusSoftCollinearLimit(const BornPhSpPt& born, double muR2) const;

        double AgluonLeg1Plus(const RealPhSpPt& born, double muR2) const;
        double AgluonLeg2Minus(const RealPhSpPt& born, double muR2) const;

        double AgluonLeg1PlusCollinearLimit(const BornPhSpPt& born, const RadiationVariables& rad, double muR2) const;
        double AgluonLeg2MinusCollinearLimit(const BornPhSpPt& born, const RadiationVariables& rad, double muR2) const;

        double AgluonLeg1PlusLimitAware(const BornPhSpPt& born, const RadiationVariables& rad, double muR2) const;
        double AgluonLeg2MinusLimitAware(const BornPhSpPt& born, const RadiationVariables& rad, double muR2) const;

        double RHatgluonLeg1Plus(const BornPhSpPt& born, const std::array<double, 3>& unitCube, double muF2, double muR2) const;
        double RHatgluonLeg2Minus(const BornPhSpPt& born, const std::array<double, 3>& unitCube, double muF2, double muR2) const;

        RadiationVariables sampleUniformRad(const BornPhSpPt& born, const std::array<double, 3>& unitCube) const;
        double unitCubeJacobian(const RadiationVariables& rad) const;

        double collinearRemnantPlusQQ(const BornPhSpPt& born, double unitX, double muF2) const;
        double collinearRemnantMinusQQ(const BornPhSpPt& born, double unitX, double muF2) const;

        double collinearRemnantGluonLeg1(const BornPhSpPt& born, double unitX, double muF2) const;
        double collinearRemnantGluonLeg2(const BornPhSpPt& born, double unitX, double muF2) const;

        double collinearRemnantContribution(const BornPhSpPt& born, const std::array<double, 4>& unitX, double muF2, double muR2) const;

        double realEndpointRemnantQQbar(const BornPhSpPt& born, const std::array<double, 3>& unitCube, double muF2, double muR2) const;

        BBarIntegrationPoint generateIntegrationPoint(
            const BornPhSpPt& born, 
            const BornChannel& bornChannel, 
            const double muF2,
            const double muR2,
            const std::array<double, 3>& unitCube
        ) const;

        double bTilde(
            const BornPhSpPt& born,
            const BornChannel& bornChannel, 
            const std::array<double, 3>& unitCube
        ) const;

    private:
        const Process& m_process;
        const Config& m_config;
        std::shared_ptr<ISRRealPhaseSpace> m_realPhaseSpace;
        
        const BornVirtualHandler m_bornVirtual;
        const CollinearRemnantsHandler m_collRemn;

        int m_sumSigns;
        int m_nEventTrials;
        double m_maxWeight;
    };

} // namespace powheg_dy
