#pragma once
 
#include "powheg_dy/base.h"
#include "powheg_dy/config/config.h"
#include "powheg_dy/math/four_vector.h"

namespace powheg_dy
{
    // TODO: Generalize
    struct BornChannel
    {
        // PDG ids of the incoming Born partons
        int id1;
        int id2;

        // TODO: Remove this
        int flavour;

        std::vector<int> outIDs;
    };

    struct BornPhSpPt
    {
        // INDEPENDENT VARIABLES

        // mass and rapidity of the intermediate gauge boson
        double mBoson = 0.0;
        double yB = 0.0;

        // angular configuration of the leptons
        double cosTh = 0.0;
        double phi = 0.0;     // everything is invariant under rotations around the beam axis. 

        BornChannel channel;

        // DEPENDENT VARIABLES

        // born-level momentum fractions and invariance mass-squared
        double x1Bar;
        double x2Bar;
        double sHat;

        // born-level momenta of the participating particles
        FourVector p1Bar;
        FourVector p2Bar;
        FourVector pBoson;
        FourVector pLMinus;
        FourVector pLPlus;

        double amp2;
        double f1;
        double f2;

        // monte carlo variables
        double jacobianOld;
        double jacobian;
        double weight;
        int sign;
    };

    class BornPhaseSpace
    {
    public:
        BornPhaseSpace(const Config& config)
            : m_config(config)
        {
        }

        BornPhSpPt samplePoint(double rands[3]) const;
        void reconstructMomenta(BornPhSpPt& point) const ;
        
    private:
        const Config& m_config;
    };


} // namespace powheg_dy
