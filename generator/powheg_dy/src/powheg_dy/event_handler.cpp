#include "event_handler.h"

#include "powheg_dy/math.h"
#include "powheg_dy/process.h"
#include "powheg_dy/born_event.h"

namespace powheg_dy
{
    namespace 
    {
        static constexpr double __ALLOWED_KINEMATIC_MISMATCH = 1.0e-10;

    } // namespace

    Event EventHandler::reconstructEvent(const PhaseSpacePoint& point, 
            const BornEvent& bornEvent, const Emission& emission) const
    {
        Event event(point, bornEvent, emission);

        // In the Born event, \cos(\theta) is the angle between the quark and the lepton. Here, 
        // we need the angle between leg 1 and the lepton. Thus, if the quark is actually on 
        // leg 2, we need to flip the sign of \cos(\theta).
        const double cosThLeg1 = bornEvent.partonId > 0 ? point.cosTh : -point.cosTh;

        const double pT = sqrt(emission.t);
        const double mT = sqrt(point.sHat + emission.t);

        const double x1PreEm = point.x1 / (emission.leg == 1 ? emission.z : 1.0);
        const double x2PreEm = point.x2 / (emission.leg == 1 ? 1.0 : emission.z);

        const double yBoson = _solveBosonRapidityFromMasslessGluon(point, emission, x1PreEm, x2PreEm, mT);

        event.pBoson = {
            mT * cosh(yBoson),
            -pT * cos(emission.phi),
            -pT * sin(emission.phi),
            mT * sinh(yBoson)
        };

        event.p1In = {
            0.5 * x1PreEm * m_process.sqrtS(),
            0.0,
            0.0,
            0.5 * x1PreEm * m_process.sqrtS()
        };

        event.p2In = {
            0.5 * x2PreEm * m_process.sqrtS(),
            0.0,
            0.0,
            -0.5 * x2PreEm * m_process.sqrtS()
        };

        const double p = sqrt(event.point.sHat) / 2.0;
        const double sinTh = sqrt(1.0 - cosThLeg1 * cosThLeg1);

        const FourVector p1Rest = {
            p,
            p * sinTh * cos(event.point.phi),
            p * sinTh * sin(event.point.phi),
            p * cosThLeg1
        };

        const FourVector p2Rest = { p1Rest.e, -p1Rest.getThreeVec() };

        event.p1Out = p1Rest.boost(event.pBoson);
        event.p2Out = p2Rest.boost(event.pBoson);
        event.pGluon = event.p1In + event.p2In - event.pBoson;

        const FourVector totalIn = event.p1In + event.p2In;
        const FourVector totalOut = event.p1Out + event.p2Out + event.pGluon;
        double mismatch = (totalIn - totalOut) * (totalIn - totalOut) / point.sHat;

        assert(abs(mismatch) < __ALLOWED_KINEMATIC_MISMATCH);

        return event;
    }

    double EventHandler::_solveBosonRapidityFromMasslessGluon(const PhaseSpacePoint& point, const Emission& emission, 
        double x1PreEm, double x2PreEm, double mT) const
        {
            if (emission.rejected)
                return point.yBoson;

            double a = x1PreEm * m_process.sqrtS();
            double b = x2PreEm * m_process.sqrtS();

            double C = (a * b + point.sHat) / mT;
            double disc = C * C - 4.0 * a * b;
            
            assert(disc >= 0);      // Unphysical kinematics

            double sqrtDisc = sqrt(disc);

            double uPlus = (C + sqrtDisc) / (2.0 * b);
            double uMinus = (C - sqrtDisc) / (2.0 * b);

            double yPlus = log(uPlus);
            double yMinus = log(uMinus);

            double yReal = log(x1PreEm / x2PreEm) / 2.0;

            // pick the solution for which the gluon moves in the same z-direction as the quark 
            // it was emitted from (in the partonic COM frame)
            if (emission.leg == 1) 
                return (yPlus < yReal) ? yPlus : yMinus;
            else
                return (yPlus > yReal) ? yPlus : yMinus;
        }

} // namespace powheg_dy
