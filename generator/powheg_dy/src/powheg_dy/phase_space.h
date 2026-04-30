#pragma once
 
#include "four_vector.h"

namespace powheg_dy
{
    class Process;

    struct BornChannel
    {
        // PDG ids of the incoming Born partons
        int id1 = 0;
        int id2 = 0;

        // absolute quark flavour: 1=d, 2=u, 3=s, 4=c, 5=b, 6=t
        int flavour = 0;
    };

    struct BornPhaseSpacePt
    {
        // INDEPENDENT VARIABLES

        // mass and rapidity of the intermediate gauge boson
        double mB = 0.0;
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
        FourVector pB;
        FourVector pLMinus;
        FourVector pLPlus;

        // monte carlo variables
        double jacobian;
        double weight;
    };

    class BornPhaseSpace
    {
    public:
        BornPhaseSpace(const Process& process)
            : m_process(process)
        {
        }

        BornPhaseSpacePt samplePoint(double rands[3]) const;
        void reconstructMomenta(BornPhaseSpacePt& point) const ;
        
    private:
        const Process& m_process;
    };


} // namespace powheg_dy
