#include "born_event_generator.h"

#include "powheg_dy/process.h"
#include "powheg_dy/math/rand.h"
#include "powheg_dy/matrix_elements.h"

namespace powheg_dy
{
    namespace
    {
        static const std::vector<int> VALID_PARTONS_ON_LEG1 = { -5, -4, -3, -2, -1, 1, 2, 3, 4, 5 };

    } // anonymous namespace

    void BornEventGenerator::computeWeightAndSampleChannel(BornPhSpPt& born) const
    {
        // For each quark flavour, compute the individual contribution to the cross section
        auto channels = computePartonChannelContributions(born);

        double totalWeight = 0.0;
        for (auto [partonId, dSigma] : channels)
            totalWeight += dSigma;

        powheg_assert(std::isfinite(totalWeight));
        powheg_assert(totalWeight >= 0.0);

        // Sample the parton channel by their relative contribution to dSigma
        double u = rand(0.0, totalWeight);
        for (const auto& [borncopy, dSigma] : channels)
        {
            if (u < dSigma)
            {   
                born = borncopy;
                born.weight = totalWeight;
                break;
            }

            u -= dSigma;
        }
    }

    std::vector<std::tuple<BornPhSpPt, double>> BornEventGenerator::computePartonChannelContributions(const BornPhSpPt& born) const
    {
        double physicsPrefactor = 1.0 / (64.0 * PI * PI * m_config.S * born.sHat);

        std::vector<std::tuple<BornPhSpPt, double>> channels;
        channels.reserve(VALID_PARTONS_ON_LEG1.size());
        
        for (int partonId : VALID_PARTONS_ON_LEG1)
        {
            BornPhSpPt borncopy = born; 
            borncopy.channel = { partonId, -partonId, abs(partonId) };
            m_bornPhaseSpace->reconstructMomenta(borncopy);

            // Compute the luminosity factors
            double f  = m_config.PDF->xfxQ2(borncopy.channel.id1, borncopy.x1Bar, borncopy.sHat) / borncopy.x1Bar;
            double fb = m_config.PDF->xfxQ2(borncopy.channel.id2, borncopy.x2Bar, borncopy.sHat) / borncopy.x2Bar;

            // Compute the event weight
            const double amp2 = m_process.bornContribution(borncopy);
            double weight = f * fb * amp2;

            channels.push_back({ borncopy, borncopy.jacobianOld * physicsPrefactor * weight });
        }

        return channels;
    }

} // namespace powheg_dy
