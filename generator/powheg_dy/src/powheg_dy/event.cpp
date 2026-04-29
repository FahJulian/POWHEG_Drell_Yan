#include "event.h"

#include "powheg_dy/process.h"
#include "powheg_dy/born_event.h"

#include <math.h>

namespace powheg_dy
{
    namespace 
    {
        static constexpr double __ALLOWED_KINEMATIC_MISMATCH = 1.0e-10;

    } // namespace

    void Event::reconstructMomenta()
    {
        // In the Born event, \cos(\theta) is the angle between the quark and the lepton. Here, 
        // we need the angle between leg 1 and the lepton. Thus, if the quark is actually on 
        // leg 2, we need to flip the sign of \cos(\theta).
        double cosThLeg1 = m_bornEvent.getPartonId() > 0 ? m_bornEvent.getCosTh() : -m_bornEvent.getCosTh();

        double pT = std::sqrt(m_emission.getT());
        double mT = std::sqrt(m_bornEvent.getS() + m_emission.getT());

        double x1PreEm = m_bornEvent.getX1() / (m_emission.getLeg() == 1 ? m_emission.getZ() : 1.0);
        double x2PreEm = m_bornEvent.getX2() / (m_emission.getLeg() == 1 ? 1.0 : m_emission.getZ());

        double yBoson = _solveBosonRapidityFromMasslessGluon(x1PreEm, x2PreEm, mT);

        m_pBoson = {
            mT * std::cosh(yBoson),
            -pT * std::cos(m_emission.getPhi()),
            -pT * std::sin(m_emission.getPhi()),
            mT * std::sinh(yBoson)
        };

        m_p1In = {
            0.5 * x1PreEm * m_process.getSqrtS(),
            0.0,
            0.0,
            0.5 * x1PreEm * m_process.getSqrtS()
        };

        m_p2In = {
            0.5 * x2PreEm * m_process.getSqrtS(),
            0.0,
            0.0,
            -0.5 * x2PreEm * m_process.getSqrtS()
        };

        double p = std::sqrt(m_bornEvent.getS()) / 2.0;
        double sinTh = std::sqrt(1.0 - cosThLeg1 * cosThLeg1);

        FourVector p1Rest = {
            p,
            p * sinTh * std::cos(m_bornEvent.getPhi()),
            p * sinTh * std::sin(m_bornEvent.getPhi()),
            p * cosThLeg1
        };

        FourVector p2Rest = { p1Rest.e, -p1Rest.getThreeVec() };

        m_p1Out = p1Rest.boost(m_pBoson);
        m_p2Out = p2Rest.boost(m_pBoson);
        m_pGluon = m_p1In + m_p2In - m_pBoson;

        ASSERT(_hasValidKinematics(), "Invalid kinematics");
    }

    std::string Event::toString() const
    {
        FourVector sum = m_p1Out + m_p2Out;
        double pT = std::sqrt(sum.pX*sum.pX + sum.pY*sum.pY);

        return std::to_string(m_bornEvent.getMBoson()) + ", "
            + std::to_string(m_bornEvent.getCosTh()) + ", "
            + std::to_string(pT);
    }

    double Event::_solveBosonRapidityFromMasslessGluon(double x1PreEm, double x2PreEm, double mT) const
    {
        if (m_emission.isRejected())
            return m_bornEvent.getYBoson();

        double a = x1PreEm * m_process.getSqrtS();
        double b = x2PreEm * m_process.getSqrtS();

        double C = (a * b + m_bornEvent.getS()) / mT;
        double disc = C * C - 4.0 * a * b;
        
        ASSERT(disc >= 0, "Unphysical kinematics");

        double sqrtDisc = std::sqrt(disc);

        double uPlus = (C + sqrtDisc) / (2.0 * b);
        double uMinus = (C - sqrtDisc) / (2.0 * b);

        double yPlus = std::log(uPlus);
        double yMinus = std::log(uMinus);

        double yReal = std::log(x1PreEm / x2PreEm) / 2.0;

        // pick the solution for which the gluon moves in the same z-direction as the quark 
        // it was emitted from (in the partonic COM frame)
        if (m_emission.getLeg() == 1) 
            return (yPlus < yReal) ? yPlus : yMinus;
        else
            return (yPlus > yReal) ? yPlus : yMinus;
    }

    bool Event::_hasValidKinematics() const
    {
        FourVector totalIn = m_p1In + m_p2In;
        FourVector totalOut = m_p1Out + m_p2Out + m_pGluon;
        double mismatch = (totalIn - totalOut) * (totalIn - totalOut) / m_bornEvent.getS();

        return std::abs(mismatch) < __ALLOWED_KINEMATIC_MISMATCH;
    }

} // namespace powheg_dy
