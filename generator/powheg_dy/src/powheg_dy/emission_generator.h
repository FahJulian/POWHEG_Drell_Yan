#pragma once

#include "powheg_dy/event.h"
#include "powheg_dy/emission.h"
#include "powheg_dy/phase_space.h"
#include "powheg_dy/real_phase_space.h"

namespace powheg_dy
{
    class Process;

    class EmissionGenerator
    {
    public:
        EmissionGenerator(const Process& process)
            : m_process(process),
              m_realPhaseSpace(process)
        {
        }

        Emission generateEmission(const BornPhSpPt& born) const;

    private:
        Emission _generateCandidate(
            const BornPhSpPt& born,
            RadiationRegion region
        ) const;

        double _exactRadiationDensity(const RealPhSpPt& born) const;
        double _upperRadiationDensity(const RealPhSpPt& born, double kt2Trial) const;

        RadiationVariables _sampleTrialRadiation(
            const BornPhSpPt& born,
            RadiationRegion region,
            double kt2Trial
        ) const;

        double _sampleTrialKt2(const BornPhSpPt& born, double ptMax2) const;
        double _sampleTrialXi(const BornPhSpPt& born, double pT2) const;
        double _computeYForRegion(const BornPhSpPt& born, double pT2, double xi, RadiationRegion region) const;
        double _sampleTrialPhi() const;

        double _VExact(double pt2, double m2, int nF) const;
        double _VTildeLog(double logPt2, double sHat, int nF) const;
        double _integrateVTildeLog(double logPt2, double logKtmax2, double mB2, int nF) const;

        double _globalKt2Max(const BornPhSpPt& born) const;

    private:
        const Process& m_process;
        FKSRealPhaseSpace m_realPhaseSpace;
    };

} // namespace powheg_dy
