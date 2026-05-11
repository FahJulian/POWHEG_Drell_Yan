#include "bbar_integrator.h"

#include "powheg_dy/process.h"
#include "powheg_dy/math/rand.h"

namespace powheg_dy
{
namespace 
{
    static constexpr double GEV_SQ_TO_PB = 0.389379338e9;
    static constexpr double SECURITY_FACTOR = 1.1;
    static constexpr int MAX_TRIALS = 10000;

} // anonymous namespace 
    
    void BBarIntegrator::computeWeightAndSampleChannel(BornPhSpPt& born) const
    {
        std::vector<BornChannel> channels = m_process.bornChannels();
        
        std::vector<std::pair<BornPhSpPt, double>> weights = {};
        weights.reserve(channels.size());

        double totalWeight = 0.0;

        for (const auto& channel : channels)
        {
            BornPhSpPt bornCopy = born;
            bornCopy.channel = channel;
            m_bornPhaseSpace->reconstructMomenta(bornCopy);

            const double weight = bTilde(bornCopy, { });

            totalWeight += weight;
            weights.push_back({ bornCopy, weight });
        }

        born.weight = totalWeight;

        // Sample the parton channel by their relative contribution to dSigma
        double u = rand(0.0, totalWeight);
        for (auto& [bornCopy, weight] : weights)
        {
            if (u < weight)
            {   
                bornCopy.weight = totalWeight;
                born = bornCopy;

                break;
            }

            u -= weight;
        }
    }

    BornPhSpPt BBarIntegrator::sampleAccordingtoBTilde()
    {
        for (int trials = 0; trials < MAX_TRIALS; trials++)
        {
            double rands[3] = { rand(), rand(), rand() };
            BornPhSpPt born = m_bornPhaseSpace->samplePoint(rands);
            
            computeWeightAndSampleChannel(born);
            assert(born.weight <= m_maxWeight, "Born weight " << born.weight << " exceeds max weight " << m_maxWeight);

            double u = rand();
            if (u < born.weight / m_maxWeight)
            {
                m_nEventTrials += trials;
                return born;
            }
        }

        assert(false, "Bbar sampling stuck in loop.");
        return { };
    }

    double BBarIntegrator::bTilde(const BornPhSpPt& born, const RadiationVariables& rad) const
    {
        const double physicsPrefactor = 1.0 / (64.0 * PI * PI * m_config.S * born.sHat);

        const double f  = m_config.PDF->xfxQ2(born.channel.id1, born.x1Bar, born.sHat) / born.x1Bar;
        const double fb = m_config.PDF->xfxQ2(born.channel.id2, born.x2Bar, born.sHat) / born.x2Bar;

        const double ampBorn = m_process.born(born);

        return physicsPrefactor * born.jacobian * f * fb * ampBorn;
    }

    void BBarIntegrator::determineMaxWeight()
    {
        Log::info("Determining Born Veto weight");

        double max_dSigma = 0.0;
        
        for (size_t i = 0; i < m_config.N_TRIAL_EVENTS; i++)
        {   
            double rands[3] = { rand(), rand(), rand() };
            BornPhSpPt born = m_bornPhaseSpace->samplePoint(rands);
            computeWeightAndSampleChannel(born);

            if (born.weight > max_dSigma)
                max_dSigma = born.weight;
        }

        m_maxWeight = SECURITY_FACTOR * max_dSigma;

        Log::info << "Done determining Born Veto weight." << std::endl << std::endl;
    }

    double BBarIntegrator::getTotalCrossSection() const 
    {
        return m_maxWeight * GEV_SQ_TO_PB * getAcceptanceRatio();
    }

    double BBarIntegrator::getAcceptanceRatio() const 
    {
        return static_cast<double>(m_config.N_ACCEPTED_EVENTS) / m_nEventTrials;
    }

    void BBarIntegrator::clear()
    {
        m_maxWeight = 0.0;
        m_nEventTrials = 0;
    }

} // namespace powheg_dy
