#include "born_event.h"

#include "powheg_dy/rand.h"
#include "powheg_dy/math.h"
#include "powheg_dy/process.h"

namespace powheg_dy
{
    namespace
    {
        static const std::vector<int> __VALID_PARTONS_ON_LEG1 = { -5, -4, -3, -2, -1, 1, 2, 3, 4, 5 };

    } // namespace

    BornEvent BornEventGenerator::computeWeightAndSampleParton(const PhaseSpacePoint& point) const
    {
        BornEvent event;

        // For each quark flavour, compute the individual contribution to the cross section
        auto channels = _computePartonChannelContributions(point);

        event.dSigma = 0.0;
        for (auto [partonId, dSigma] : channels)
            event.dSigma += dSigma;

        // Sample the parton channel by their relative contribution to dSigma
        double u = rand(0.0, event.dSigma);
        for (auto [partonId, dSigma] : channels)
        {
            event.partonId = partonId;
            if (u < dSigma)
                break;

            u -= dSigma;
        }

        return event;
    }

    std::vector<std::tuple<int, double>> BornEventGenerator::_computePartonChannelContributions(const PhaseSpacePoint& point) const
    {
        double physicsPrefactor = m_process.ALPHA() * m_process.ALPHA() / 2.0 / m_process.NC() 
            / m_process.sqrtS() / m_process.sqrtS() / point.mBoson;

        std::vector<std::tuple<int, double>> channels;
        channels.reserve(__VALID_PARTONS_ON_LEG1.size());
        
        for (int partonId : __VALID_PARTONS_ON_LEG1)
        {
            // Compute the luminosity factors
            double f  = m_process.getPdfs()->xfxQ2( partonId, point.x1, point.sHat) / point.x1;
            double fb = m_process.getPdfs()->xfxQ2(-partonId, point.x2, point.sHat) / point.x2;

            // Compute the coupling factors
            bool upType = std::abs(partonId) % 2 == 0;
            auto [c1, c2] = _neutralCurrentCouplingFactors(upType, point.sHat);

            // Compute the event weight
            double weight = f * fb * (c1 * (1.0 + point.cosTh * point.cosTh) + 2.0 * c2 * point.cosTh);

            channels.push_back({ partonId, point.invSamplingFact * physicsPrefactor * weight });
        }

        return channels;
    }

    std::tuple<double, double> BornEventGenerator::_neutralCurrentCouplingFactors(bool upType, double mSq) const
        {
            // quark charge and axial and vector couplings
            double qQ = upType ? 2.0 / 3.0 : -1.0 / 3.0;
            double aQ = upType ? 0.5 : -0.5;
            double vQ = upType ? 0.5 - (4.0 / 3.0) * m_process.S_W_SQ() : -0.5 + (2.0 / 3.0) * m_process.S_W_SQ();

            double dZ = (mSq - m_process.M_Z() * m_process.M_Z()) * (mSq - m_process.M_Z() * m_process.M_Z()) 
                + m_process.M_Z() * m_process.M_Z() * m_process.GAMMA_Z() * m_process.GAMMA_Z();

            double ReChi   = m_process.KAPPA() * mSq * (mSq - m_process.M_Z() * m_process.M_Z()) / dZ;
            double AbsChiSq = m_process.KAPPA() * m_process.KAPPA() * mSq * mSq / dZ;

            double hU = qQ * qQ
                    - 2.0 * qQ * m_process.V_L() * vQ * ReChi
                    + (m_process.V_L() * m_process.V_L() + m_process.A_L() * m_process.A_L()) * (vQ * vQ + aQ * aQ) * AbsChiSq;

            double hF = -2.0 * qQ * m_process.A_L() * aQ * ReChi
                    + 4.0 * m_process.V_L() * m_process.A_L() * vQ * aQ * AbsChiSq;

            return {hU, hF};
        }

} // namespace powheg_dy
