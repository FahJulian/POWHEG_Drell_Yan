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
        BornPhSpPt sampleAccordingtoBTilde();

        double getAcceptanceRatio() const;
        double getTotalCrossSection() const;

    private:
        struct Counterterms
        {
            double qqbar;
            double gluonLeg1;
            double gluonLeg2;
        };

        bool computeWeightAndSampleChannel(BornPhSpPt& born) const;

        double bTilde(const BornPhSpPt& born, const std::array<double, 3>& unitCube) const;
        Counterterms counterterms(const BornPhSpPt& born, const RadiationVariables& rad) const;

        double bornContribution(const BornPhSpPt& born, double muF2) const;
        double bornPlusVirtualContribution(const BornPhSpPt& born) const;
        double realMinusCounterTermContribution(const BornPhSpPt& born, const std::array<double, 3>& unitCube) const;

        double luminosity(const BornPhSpPt& born, const RadiationVariables& rad, int id1, int id2, double muF2) const;

        double RHatqqbarPlus(const BornPhSpPt& born, const std::array<double, 3>& unitCube, double muF2, double muR2) const;
        double RHatqqbarMinus(const BornPhSpPt& born, const std::array<double, 3>& unitCube, double muF2, double muR2) const;

        double AqqbarPlusLimitAware(const BornPhSpPt& born, const RadiationVariables& rad, double muF2, double muR2) const;
        double AqqbarMinusLimitAware(const BornPhSpPt& born, const RadiationVariables& rad, double muF2, double muR2) const;
        
        double AqqbarPlus(const RealPhSpPt& real, double muF2, double muR2) const;
        double AqqbarMinus(const RealPhSpPt& real, double muF2, double muR2) const;

        double AqqbarPlusCollinearLimit(const BornPhSpPt& born, const RadiationVariables& rad, double muF2, double muR2) const;
        double AqqbarMinusCollinearLimit(const BornPhSpPt& born, const RadiationVariables& rad, double muF2, double muR2) const;

        double AqqbarPlusSoftLimit(const BornPhSpPt& born, const RadiationVariables& rad, double muF2, double muR2) const;
        double AqqbarMinusSoftLimit(const BornPhSpPt& born, const RadiationVariables& rad, double muF2, double muR2) const;

        double AqqbarPlusSoftCollinearLimit(const BornPhSpPt& born, const RadiationVariables& rad, double muF2, double muR2) const;
        double AqqbarMinusSoftCollinearLimit(const BornPhSpPt& born, const RadiationVariables& rad, double muF2, double muR2) const;

        double AgluonLeg1Plus(const RealPhSpPt& born, double muF2, double muR2) const;
        double AgluonLeg2Minus(const RealPhSpPt& born, double muF2, double muR2) const;

        double AgluonLeg1PlusCollinearLimit(const BornPhSpPt& born, const RadiationVariables& rad, double muF2, double muR2) const;
        double AgluonLeg2MinusCollinearLimit(const BornPhSpPt& born, const RadiationVariables& rad, double muF2, double muR2) const;

        double AgluonLeg1PlusLimitAware(const BornPhSpPt& born, const RadiationVariables& rad, double muF2, double muR2) const;
        double AgluonLeg2MinusLimitAware(const BornPhSpPt& born, const RadiationVariables& rad, double muF2, double muR2) const;

        double RHatgluonLeg1Plus(const BornPhSpPt& born, const std::array<double, 3>& unitCube, double muF2, double muR2) const;
        double RHatgluonLeg2Minus(const BornPhSpPt& born, const std::array<double, 3>& unitCube, double muF2, double muR2) const;

        RadiationVariables sampleUniformRad(const BornPhSpPt& born, const std::array<double, 3>& unitCube) const;
        double unitCubeJacobian(const RadiationVariables& rad) const;

    private:
        const Process& m_process;
        const Config& m_config;
        std::shared_ptr<BornPhaseSpace> m_bornPhaseSpace;
        std::shared_ptr<FKSRealPhaseSpace> m_realPhaseSpace;

        int m_nEventTrials;
        double m_maxWeight;
    };

} // namespace powheg_dy
