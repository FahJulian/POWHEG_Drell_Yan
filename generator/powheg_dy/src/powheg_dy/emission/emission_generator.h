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

        Emission generateEmission(const BornPhSpPt& born) const;

    private:
        double _upperRadiationDensity(const RealPhSpPt& born, double kt2Trial) const;

        RadiationVariables _sampleTrialRadiation(const BornPhSpPt& born, double kt2Trial) const;

        double _sampleTrialKt2(const BornPhSpPt& born, double ptMax2) const;
        double _sampleTrialXi(const BornPhSpPt& born, double pT2) const;
        double _sampleY(const BornPhSpPt& born, double pT2, double xi) const;
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
