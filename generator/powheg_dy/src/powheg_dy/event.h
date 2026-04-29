#pragma once

#include "powheg_dy/base.h"
#include "powheg_dy/emission.h"
#include "powheg_dy/born_event.h"
#include "powheg_dy/four_vector.h"

#include <string>

namespace powheg_dy
{
    class Process;

    class Event
    {
    public:
        Event(const Process& process, const BornEvent& bornEvent, const Emission& emission)
            : m_process(process), m_bornEvent(bornEvent), m_emission(emission)
        {
        }

        void reconstructMomenta();
        std::string toString() const;
        const BornEvent& getBornEvent() const { return m_bornEvent; }
        const Emission& getEmission() const { return m_emission; }
        const FourVector& getP1In() const { return m_p1In; }
        const FourVector& getP2In() const { return m_p2In; }
        const FourVector& getP1Out() const { return m_p1Out; }
        const FourVector& getP2Out() const { return m_p2Out; }
        const FourVector& getPBoson() const { return m_pBoson; }
        const FourVector& getPGluon() const { return m_pGluon; }

    private:
        double _solveBosonRapidityFromMasslessGluon(double x1PreEm, double x2PreEm, double mT) const;
        bool _hasValidKinematics() const;

    private:
        const Process& m_process;
        const BornEvent m_bornEvent;
        Emission m_emission;
        FourVector m_p1In;
        FourVector m_p2In;
        FourVector m_p1Out;
        FourVector m_p2Out;
        FourVector m_pBoson;
        FourVector m_pGluon;
    };
    
} // namespace powheg_dy
