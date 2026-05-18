#pragma once

#include "powheg/base.h"
#include "powheg/config/config.h"
#include "powheg/phase_space/born_phase_space_point.h"
#include "powheg/phase_space/real_phase_space.h"

namespace powheg
{
    class EmissionSampler
    {
    public:
        EmissionSampler(const Config& config, std::shared_ptr<ISRRealPhaseSpace> realPhaseSpace)
            : m_config(config),
              m_realPhaseSpace(std::move(realPhaseSpace))
        {
        }

        RadiationVariables sampleTrialRadiation(const BornPhSpPt& born, double kt2Trial) const;
        double globalKt2Max(const BornPhSpPt& born) const;
        double sampleTrialKt2(const BornPhSpPt& born, double ptMax2, double& logR) const;
        double upperRadiationDensity(const RadiationVariables& rad, double kt2Trial) const;

    private:
        double sampleTrialXi(const BornPhSpPt& born, double pT2) const;
        double sampleY(const BornPhSpPt& born, double pT2, double xi) const;
        double sampleTrialPhi() const;
        
        double integrateVTilde(double pt2, double kt2max, double sBorn, double lambda2, double nQ, double beta0) const;
        double vExactOverVTilde(const BornPhSpPt& born, double pt2) const;


    private:
        const Config& m_config;
        std::shared_ptr<ISRRealPhaseSpace> m_realPhaseSpace;
    };

} // namespace powheg
