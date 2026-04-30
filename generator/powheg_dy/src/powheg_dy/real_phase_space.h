#pragma once

#include "powheg_dy/phase_space.h"
#include "powheg_dy/four_vector.h"

namespace powheg_dy
{
    class Process;

    struct RadiationVariables
    {
        // FKS/POWHEG radiation variables
        double xi  = 0.0;   // emitted-parton energy fraction in the real partonic CM
        double y   = 0.0;   // cos(theta) relative to leg-1 direction
        double phi = 0.0;
    };

    struct RealPhaseSpacePt
    {
        BornPhaseSpacePt underlyingBorn;
        RadiationVariables rad;

        double x1 = 0.0;
        double x2 = 0.0;
        double sHatReal = 0.0;

        // POWHEG/FKS radiation quantities
        double kt2 = 0.0;
        double radJacobian = 0.0;

        FourVector p1In;
        FourVector p2In;

        FourVector pLMinus;
        FourVector pLPlus;
        FourVector pBoson;
        FourVector pRadiated;
    };

    class FKSRealPhaseSpace
    {
    public:
        FKSRealPhaseSpace(const Process& process)
            : m_process(process)
        {
        }

        RealPhaseSpacePt reconstruct(
            const BornPhaseSpacePt& born, 
            const RadiationVariables& rad
        ) const;

        double xiMax(
            const BornPhaseSpacePt& born,
            double y
        ) const;

    private:
        const Process& m_process;
    };

} // namespace powheg_dy