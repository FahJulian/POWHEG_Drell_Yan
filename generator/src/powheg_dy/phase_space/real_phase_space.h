#pragma once

#include "powheg_dy/base.h"
#include "powheg_dy/config.h"
#include "powheg_dy/math/four_vector.h"
#include "powheg_dy/phase_space/born_phase_space.h"

namespace powheg_dy
{
    struct RadiationVariables
    {
        // FKS/POWHEG radiation variables
        double xi  = 0.0;   // emitted-parton energy fraction in the real partonic CM
        double y   = 0.0;   // cos(theta) relative to leg-1 direction
        double phi = 0.0;
        double xiMax = 0.0;
    };

    struct RealChannel
    {
        // PDG ids of the incoming legs
        int id1;
        int id2;

        // PDG ids of the outgoing legs
        std::vector<int> outIDs;
    };

    struct RealPhSpPt
    {
        BornPhSpPt underlyingBorn;
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
        FKSRealPhaseSpace(const Config& config)
            : m_config(config)
        {
        }

        RealPhSpPt reconstruct(
            const BornPhSpPt& born, 
            const RadiationVariables& rad
        ) const;

        double xiMax(
            const BornPhSpPt& born,
            double y
        ) const;

    private:
        const Config& m_config;
    };

} // namespace powheg_dy