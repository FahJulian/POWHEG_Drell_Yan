#include "born_event.h"

#include "powheg_dy/rand.h"
#include "powheg_dy/math.h"
#include "powheg_dy/process.h"
#include "powheg_dy/matrix_elements.h"

#include <cmath>

namespace powheg_dy
{
    namespace
    {
        static const std::vector<int> __VALID_PARTONS_ON_LEG1 = { -5, -4, -3, -2, -1, 1, 2, 3, 4, 5 };

    } // namespace

    void BornEventGenerator::computeWeightAndSampleChannel(BornPhaseSpacePt& point) const
    {
        // For each quark flavour, compute the individual contribution to the cross section
        auto channels = _computePartonChannelContributions(point);

        point.weight = 0.0;
        for (auto [partonId, dSigma] : channels)
            point.weight += dSigma;

        assert(point.weight >= 0.0 && std::isfinite(point.weight));

        // Sample the parton channel by their relative contribution to dSigma
        double u = rand(0.0, point.weight);
        for (auto [channel, dSigma] : channels)
        {
            point.channel = channel;
            if (u < dSigma)
                break;

            u -= dSigma;
        }
    }

    std::vector<std::tuple<BornChannel, double>> BornEventGenerator::_computePartonChannelContributions(const BornPhaseSpacePt& point) const
    {
        double physicsPrefactor = m_process.ALPHA() * m_process.ALPHA() / 2.0 / m_process.NC() 
            / m_process.sqrtS() / m_process.sqrtS() / point.mB;

        std::vector<std::tuple<BornChannel, double>> channels;
        channels.reserve(__VALID_PARTONS_ON_LEG1.size());
        
        for (int partonId : __VALID_PARTONS_ON_LEG1)
        {
            BornChannel channel = { partonId, -partonId, abs(partonId) };

            // Compute the luminosity factors
            double f  = m_process.getPdfs()->xfxQ2(channel.id1, point.x1Bar, point.sHat) / point.x1Bar;
            double fb = m_process.getPdfs()->xfxQ2(channel.id2, point.x2Bar, point.sHat) / point.x2Bar;

            // Compute the event weight
            double weight = f * fb * MatrixElements::bornAngularFactor(m_process, channel.flavour, point.sHat, point.cosTh);

            channels.push_back({ channel, point.jacobian * physicsPrefactor * weight });
        }

        return channels;
    }

} // namespace powheg_dy
