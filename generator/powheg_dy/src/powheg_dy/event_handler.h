#pragma once
 
#include "powheg_dy/emission.h"
#include "powheg_dy/born_event.h"
#include "powheg_dy/phase_space.h"
#include "powheg_dy/four_vector.h"

#include <string>

namespace powheg_dy
{
    class Process;

    struct Event
    {
        Event(const BornPhaseSpacePt& point, const Emission& emission)
            : point(point), emission(emission)
        {
        }

        const BornPhaseSpacePt point;
        const Emission emission;
        FourVector p1In;
        FourVector p2In;
        FourVector p1Out;
        FourVector p2Out;
        FourVector pBoson;
        FourVector pGluon;
    };

    class EventHandler
    {
    public:
        EventHandler(const Process& process)
            : m_process(process)
        {
        }

        Event reconstructEvent(const BornPhaseSpacePt& point, const Emission& emission) const;

    private:
        double _solveBosonRapidityFromMasslessGluon(const BornPhaseSpacePt& point, const Emission& emission, double x1PreEm, double x2PreEm, double mT) const;

    private:
        const Process& m_process;
    };
    
} // namespace powheg_dy
