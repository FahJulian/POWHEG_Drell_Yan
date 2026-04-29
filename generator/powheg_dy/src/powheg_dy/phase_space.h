#pragma once
 

namespace powheg_dy
{
    class Process;

    struct PhaseSpacePoint
    {
        double mBoson;
        double yBoson;
        double cosTh;

        // everything is invariant under rotations around the beam axis. 
        const double phi = 0.0;

        double x1;
        double x2;
        double sHat;

        double invSamplingFact;
    };

    class PhaseSpaceSampler
    {
    public:
        PhaseSpaceSampler(const Process& process)
            : m_process(process)
        {
        }

        PhaseSpacePoint samplePoint(double rands[3]) const;
        
    private:
        const Process& m_process;
    };


} // namespace powheg_dy
