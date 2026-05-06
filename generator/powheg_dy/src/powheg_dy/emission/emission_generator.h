#pragma once

#include "powheg_dy/event.h"
#include "powheg_dy/emission/emission.h"
#include "powheg_dy/phase_space/born_phase_space.h"
#include "powheg_dy/phase_space/real_phase_space.h"
#include "powheg_dy/matrix_elements/matrix_elements.h"

#include <vector>

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

        Emission generateEmission(const BornPhSpPt& born, int region) const;
        
    private:
        Emission _generateEmissionISR(const BornPhSpPt& born) const;

        double _upperRadiationDensity(const RealPhSpPt& born, double kt2Trial) const;

        RadiationVariables _sampleTrialRadiation(const BornPhSpPt& born, double kt2Trial) const;

        double _sampleTrialKt2(const BornPhSpPt& born, double ptMax2, double& logR) const;
        double _sampleTrialXi(const BornPhSpPt& born, double pT2) const;
        double _sampleY(const BornPhSpPt& born, double pT2, double xi) const;
        double _sampleTrialPhi() const;

        double _VExact(double pt2, double m2, int nF) const;
        double _VTildeLog(const BornPhSpPt& born, double logPt2, int nF) const;
        double _integrateVTilde(
            double pt2,
            double kt2max,
            double sBorn,
            double lambda2,
            double nQ,
            double beta0
        ) const;
        double _VExactOverVTilde(const BornPhSpPt& born, double pt2) const;

        double _globalKt2Max(const BornPhSpPt& born) const;

    private:
        const Process& m_process;
        FKSRealPhaseSpace m_realPhaseSpace;
    };

} // namespace powheg_dy
