# pragma once

#include "powheg_dy/base.h"

namespace powheg_dy
{
    class Process;
    class BornEvent;
    
    class Emission
    {
    private:
        Emission() = delete;
        Emission(const Process& process, const BornEvent& bornEvent, int leg)
            : m_process(process), m_bornEvent(bornEvent), m_leg(leg)
        {
        }

    public:
        Emission reject() { m_t = 0.0; m_z = 1.0; m_phi = 0.0; m_rejected = true; return *this; }
        bool isRejected() const { return m_rejected; }

        int getLeg() const { return m_leg; }
        double getT() const { return m_t; }
        double getZ() const { return m_z; }
        double getPhi() const { return m_phi; }

    private:
        const Process& m_process;
        const BornEvent& m_bornEvent;
        const int m_leg;
        double m_t = 0.0;
        double m_z = 0.0;
        double m_phi = 0.0;
        bool m_rejected = false;
    
    public:    
        static Emission generateFirstEmission(const Process& process, const BornEvent& bornEvent);

    private:
        static Emission _generateEmissionOnLeg(const Process& process, const BornEvent& bornEvent, int leg);
    };

} // namespace powheg_dy
