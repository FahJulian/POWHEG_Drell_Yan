# pragma once
 
#include "powheg_dy/born_event.h"
#include "powheg_dy/phase_space.h"

namespace powheg_dy
{
    class Process;
    class BornEventGenerator;
    
    struct Emission
    {
        int leg;
        double t;
        double z;
        double phi;

        bool rejected = false;

        Emission reject() { t = 0.0; z = 1.0; phi = 0.0; rejected = true; return *this; }
    };

    class EmissionGenerator
    {
    public:
        EmissionGenerator(const Process& process)
            : m_process(process)
        {
        }

    public:
        Emission generateEmission(const PhaseSpacePoint& point, const BornEvent& bornEvent) const;

    private:
        const Process& m_process;
    
    private:
        Emission _generateEmissionOnLeg(const PhaseSpacePoint& point, const BornEvent& bornEvent, int leg) const;
    };

} // namespace powheg_dy
