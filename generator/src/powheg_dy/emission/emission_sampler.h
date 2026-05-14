#pragma once

#include "powheg_dy/base.h"
#include "powheg_dy/config/config.h"
#include "powheg_dy/phase_space/born_phase_space.h"
#include "powheg_dy/phase_space/real_phase_space.h"

namespace powheg_dy
{
    class EmissionSampler
    {
    public:
        EmissionSampler(const Config& config, std::shared_ptr<FKSRealPhaseSpace> realPhaseSpace)
            : m_config(config),
              m_realPhaseSpace(std::move(realPhaseSpace))
        {
        }

        RadiationVariables sampleTrialRadiation(const BornPhSpPt& born, double kt2Trial) const;
        double globalKt2Max(const BornPhSpPt& born) const;
        double sampleTrialKt2(const BornPhSpPt& born, double ptMax2, double& logR) const;
        double upperRadiationDensity(const RealPhSpPt& born, double kt2Trial) const;

    private:
        double sampleTrialXi(const BornPhSpPt& born, double pT2) const;
        double sampleY(const BornPhSpPt& born, double pT2, double xi) const;
        double sampleTrialPhi() const;
        
        double integrateVTilde(double pt2, double kt2max, double sBorn, double lambda2, double nQ, double beta0) const;
        double vExactOverVTilde(const BornPhSpPt& born, double pt2) const;


    private:
        const Config& m_config;
        std::shared_ptr<FKSRealPhaseSpace> m_realPhaseSpace;
    };

} // namespace powheg_dy
