#include "emission_generator.h"

#include "powheg/process.h"
#include "powheg/math/rand.h"

namespace powheg
{
namespace 
{
    static constexpr double EPS_XI     = 1.0e-12;
    static constexpr int MAX_TRIALS = 100000;

    RealChannel chooseChannel(const RealOverBornContributions& contributions)
    {
        double u = rand(0.0, contributions.total);

        for (const auto& [channel, contribution] : contributions.channels)
        {
            if (u < contribution)
                return channel;
            
            u -= contribution;
        }

        throw std::runtime_error("Invalid RealOverBornContributions.");
        return { };
    }

} // anonymous namespace

    Emission EmissionGenerator::generateEmission(
        const BornPhSpPt& born,
        const BornChannel& bornChannel
    ) const
    {
        // TODO: Implement FSR with highest-bid procedure 
        
        const double amp2Born = m_process.bornAmp2(born, bornChannel);
        return generateISREmission(born, bornChannel, amp2Born);
    }

    Emission EmissionGenerator::generateISREmission(
        const BornPhSpPt& born,
        const BornChannel& bornChannel,
        const double amp2Born
    ) const
    {
        double kt2Max = m_sampler.globalKt2Max(born);

        double logR = 0;
        for (int trial = 0; trial < MAX_TRIALS && kt2Max > m_config.PT_SQ_CUTOFF; ++trial)
        {
            const double kt2Trial = m_sampler.sampleTrialKt2(born, kt2Max, logR);
            if (kt2Trial < m_config.PT_SQ_CUTOFF)
                return Emission().reject();
            
            const RadiationVariables rad = m_sampler.sampleTrialRadiation(born, kt2Trial);

            if (rad.xi <= EPS_XI || rad.xi >= 1.0 - EPS_XI || rad.xi >= rad.xiMax)
            {
                kt2Max = kt2Trial;
                continue;
            }

            const double muF2 = kt2Trial;
            const double muR2 = kt2Trial;

            const RealPhSpPt real = m_realPhaseSpace->reconstruct(born, rad);

            const RealOverBornContributions contributions 
                = getRealOverBornContributions(real, born, bornChannel, amp2Born, muF2, muR2);

            const double accRatio = contributions.total / m_sampler.upperRadiationDensity(rad, kt2Trial);

            if (accRatio > 1.0)
                Log::warn << "Acceptance ratio " << accRatio << " is greater than one, accepting emission." << Log::endl;

            if (rand() < accRatio)
            {
                RealChannel channel = chooseChannel(contributions);
                return { rad, channel, kt2Trial, false };
            }
            else
                kt2Max = kt2Trial;
        }
        
        return Emission().reject();
    }

    RealOverBornContributions EmissionGenerator::getRealOverBornContributions(
        const RealPhSpPt& real,
        const BornPhSpPt& born,
        const BornChannel& bornChannel,
        const double amp2Born,
        const double muF2,
        const double muR2
    ) const
    {
        RealOverBornContributions contributions;

        const double lumBorn = m_config.PDF->xfxQ2(bornChannel.id1, born.x1Bar, muF2) / born.x1Bar
            * m_config.PDF->xfxQ2(bornChannel.id2, born.x2Bar, muF2) / born.x2Bar;

        for (const auto& realChannel: m_process.realChannels(bornChannel))
        {
            const double lumReal = m_config.PDF->xfxQ2(realChannel.id1, real.x1, muF2) / real.x1
                * m_config.PDF->xfxQ2(realChannel.id2, real.x2, muF2) / real.x2;

            const double realAmp2 = m_process.realAmp2(real, realChannel, m_config.alphaSCMW(muR2));
            const double realOverBornPartonic = real.radJacobian * born.sHat / real.sHatReal * realAmp2 / amp2Born;
            const double realOverBorn = lumReal / lumBorn * realOverBornPartonic;

            contributions.total += realOverBorn;
            contributions.channels.push_back({ realChannel, realOverBorn });
        }

        return contributions;
    }

} // namespace powheg
