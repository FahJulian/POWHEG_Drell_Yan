#pragma once

#include "powheg/base.h"
#include "powheg/config/config.h"
#include "powheg/math/four_vector.h"
#include "powheg/phase_space/born_phase_space_point.h"

namespace powheg
{
    struct RadiationVariables
    {
        // TODO: This is used in various places, but should maybe be different structs
        double xi  = 0.0;   // emitted-parton energy fraction in the real partonic CM
        double y   = 0.0;   // cos(theta) relative to leg-1 direction
        double phi = 0.0;
        double xiMax = 0.0;
    };

    struct RealChannel
    {
        int id1;
        int id2;

        std::vector<int> outIDs;
    };

    struct RealPhSpPt
    {
        double x1 = 0.0;
        double x2 = 0.0;
        double sHatReal = 0.0;

        double radJacobian = 0.0;

        FourVector p1In;
        FourVector p2In;

        std::vector<FourVector> pOut;
    };

    class ISRRealPhaseSpace
    {
    public:
        ISRRealPhaseSpace(const Config& config)
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

} // namespace powheg