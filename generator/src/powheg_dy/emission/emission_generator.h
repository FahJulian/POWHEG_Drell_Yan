#pragma once

#include "powheg_dy/base.h"
#include "powheg_dy/config.h"
#include "powheg_dy/emission/emission.h"
#include "powheg_dy/phase_space/born_phase_space.h"
#include "powheg_dy/phase_space/real_phase_space.h"

namespace powheg_dy
{
    class Process;

    struct RealOverBornContributions
    {
        std::vector<std::pair<RealChannel, double>> channels = { };
        double total = 0.0;
    };

    class EmissionGenerator
    {
    public:
        EmissionGenerator(const Process& process, const Config& config, std::shared_ptr<FKSRealPhaseSpace> realPhaseSpace)
            : m_process(process),
              m_config(config),
              m_realPhaseSpace(std::move(realPhaseSpace))
        {
        }

        Emission generateEmission(const BornPhSpPt& born, int region) const;
        
    private:
        Emission generateISREmission(const BornPhSpPt& born) const;

        RealOverBornContributions getRealOverBornContributions(
            const RealPhSpPt& real,
            const BornPhSpPt& born,
            const double amp2Born,
            const double muF2,
            const double muR2
        ) const;
        
        RadiationVariables sampleTrialRadiation(const BornPhSpPt& born, double kt2Trial) const;
        double sampleTrialKt2(const BornPhSpPt& born, double ptMax2, double& logR) const;
        double sampleTrialXi(const BornPhSpPt& born, double pT2) const;
        double sampleY(const BornPhSpPt& born, double pT2, double xi) const;
        double sampleTrialPhi() const;
        
        double upperRadiationDensity(const RealPhSpPt& born, double kt2Trial) const;
        double integrateVTilde(double pt2, double kt2max, double sBorn, double lambda2, double nQ, double beta0) const;
        double vExactOverVTilde(const BornPhSpPt& born, double pt2) const;

        double globalKt2Max(const BornPhSpPt& born) const;

    private:
        const Process& m_process;
        const Config& m_config;
        std::shared_ptr<FKSRealPhaseSpace> m_realPhaseSpace;
    };

} // namespace powheg_dy
