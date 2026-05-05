#include "born_event_generator.h"

#include "powheg_dy/math/rand.h"
#include "powheg_dy/math/math.h"
#include "powheg_dy/process.h"
#include "powheg_dy/matrix_elements/matrix_elements.h"

#include <cmath>
#include <iostream>

namespace powheg_dy
{
    namespace
    {
        static const std::vector<int> __VALID_PARTONS_ON_LEG1 = { -5, -4, -3, -2, -1, 1, 2, 3, 4, 5 };

    } // namespace

    void BornEventGenerator::computeWeightAndSampleChannel(BornPhSpPt& point) const
    {
        // For each quark flavour, compute the individual contribution to the cross section
        auto channels = _computePartonChannelContributions(point);

        point.weight = 0.0;
        for (auto [partonId, dSigma] : channels)
            point.weight += dSigma;

        assert(std::isfinite(point.weight));
        assert(point.weight >= 0.0);

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

    std::vector<std::tuple<BornChannel, double>> BornEventGenerator::_computePartonChannelContributions(const BornPhSpPt& point) const
    {
        double physicsPrefactor = 1.0 / (64.0 * PI * PI * m_process.S() * point.sHat);

        std::vector<std::tuple<BornChannel, double>> channels;
        channels.reserve(__VALID_PARTONS_ON_LEG1.size());
        
        for (int partonId : __VALID_PARTONS_ON_LEG1)
        {
            BornPhSpPt born2 = point; 
            born2.channel = { partonId, -partonId, abs(partonId) };

            // Compute the luminosity factors
            double f  = m_process.getPdfs()->xfxQ2(born2.channel.id1, point.x1Bar, point.sHat) / point.x1Bar;
            double fb = m_process.getPdfs()->xfxQ2(born2.channel.id2, point.x2Bar, point.sHat) / point.x2Bar;

            // Compute the event weight
            const double amp2 = MatrixElements::born(m_process, born2);
            double weight = f * fb * amp2;

            channels.push_back({ born2.channel, point.jacobian * physicsPrefactor * weight });
        }

        return channels;
    }

} // namespace powheg_dy
