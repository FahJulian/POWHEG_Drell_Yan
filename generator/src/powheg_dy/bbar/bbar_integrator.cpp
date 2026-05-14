#include "bbar_integrator.h"

#include "powheg_dy/process.h"
#include "powheg_dy/math/rand.h"
#include "powheg_dy/bbar/singular_region.h"
#include "powheg_dy/bbar/collinear_remnant_channel.h"

namespace powheg_dy
{
namespace 
{
    static constexpr double GEV_SQ_TO_PB = 0.389379338e9;
    static constexpr double SECURITY_FACTOR = 2.5;
    static constexpr int MAX_TRIALS = 10000;

    struct WeightedChannel
    {
        BornChannel channel;
        double amp2Born;
        int weightSign;
        double absoluteWeight;
    };

} // anonymous namespace 
    
    BornEvent BBarIntegrator::computeWeightAndSampleChannel(const BornPhSpPt& born) const
    {
        const std::array<double, 3> unitCube = { rand(), rand(), rand() };

        std::vector<BornChannel> channels = m_process.bornChannels();
        
        std::vector<WeightedChannel> weights = {};
        weights.reserve(channels.size());
        
        double totalAbsoluteWeight = 0.0;
        for (const auto& channel : channels)
        {
            const double amp2Born = m_process.bornAmp2(born, channel);

            const double weight = bTilde(born, amp2Born, channel, unitCube);
            const double absoluteWeight = std::abs(weight);
            const int sign = weight >= 0;

            totalAbsoluteWeight += absoluteWeight;
            weights.push_back({ channel, amp2Born, sign, absoluteWeight });
        }

        // Sample the parton channel by their relative contribution to dSigma
        double u = rand(0.0, totalAbsoluteWeight);
        for (const WeightedChannel& weightedChannel : weights)
        {
            if (u < weightedChannel.absoluteWeight)
                return { born, weightedChannel.channel, weightedChannel.amp2Born, totalAbsoluteWeight, weightedChannel.weightSign };

            u -= weightedChannel.absoluteWeight;
        }

        throw std::runtime_error("Invalid Weights in BBar Sampling");
        return { };
    }

    BornEvent BBarIntegrator::sampleAccordingtoBTilde()
    {
        for (int trials = 1; trials <= MAX_TRIALS; trials++)
        {
            const BornEvent bornEvent = computeWeightAndSampleChannel(m_process.sampleBorn());

            if (bornEvent.absoluteWeight > m_maxWeight)
                Log::warn << "Born weight " << bornEvent.absoluteWeight << " exceeds max weight " << m_maxWeight << ", accepting event." << std::endl;

            const double u = rand();
            if (u < bornEvent.absoluteWeight / m_maxWeight)
            {
                m_nEventTrials += trials;
                m_sumSigns += bornEvent.weightSign;
                return bornEvent;
            }
        }

        powheg_assert(false, "Bbar sampling stuck in loop.");
        return { };
    }

    double BBarIntegrator::bTilde(
        const BornPhSpPt& born,
        const double amp2Born, 
        const BornChannel& bornChannel, 
        const std::array<double, 3>& unitCube
    ) const
    {
        const double muF2 = born.sHat;
        const double muR2 = born.sHat;

        const BBarIntegrationPoint point = generateIntegrationPoint(born, bornChannel, amp2Born, muF2, muR2, unitCube);

        double dSigma = m_bornVirtual.dSigmaBorn(point);

        if (m_config.BTILDE_BORNONLY)
            return dSigma;
        
        dSigma += m_bornVirtual.dSigmaVirtual(point, muR2);

        for (const auto& realChannel : m_process.realChannels(bornChannel))
        {
            for (const auto& singularRegion : findSingularRegions(bornChannel, realChannel))
            {
                if (const auto& collinearChannel = remnantChannelFromRegion(bornChannel, realChannel, singularRegion); 
                    collinearChannel.has_value())
                {
                    if (collinearChannel->leg == 1)
                        dSigma += m_collRemn.dSigmaCollinearRemnantsLeg1(point, collinearChannel->splitting, muF2);
                    else if (collinearChannel->leg == 2)
                        dSigma += m_collRemn.dSigmaCollinearRemnantsLeg2(point, collinearChannel->splitting, muF2);
                    else
                        throw std::runtime_error("Invalid collinear remnant leg");
                }
            }

            // TODO: Add real term
        }

        return dSigma;
    }

    void BBarIntegrator::determineMaxWeight()
    {
        
        if (m_config.BORN_VETO_WEIGHT != -1.0)
        {
            m_maxWeight = m_config.BORN_VETO_WEIGHT;
            Log::info << "Using manual Born Veto weight " << m_maxWeight << std::endl;
        }
        else
        {
            Log::info("Determining Born Veto weight");

            double maxWeight = 0.0;
            for (size_t i = 0; i < m_config.N_TRIAL_EVENTS; i++)
            {   
                const BornEvent bornEvent = computeWeightAndSampleChannel(m_process.sampleBorn());

                if (bornEvent.absoluteWeight > maxWeight)
                    maxWeight = bornEvent.absoluteWeight;
            }

            m_maxWeight = SECURITY_FACTOR * maxWeight;

            Log::info << "Done determining Born Veto weight." << std::endl << std::endl;
        }
    }

    double BBarIntegrator::getTotalCrossSection() const 
    {
        return m_maxWeight * GEV_SQ_TO_PB * static_cast<double>(m_sumSigns) / m_nEventTrials;
    }

    double BBarIntegrator::getAbsCrossSection() const 
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
        m_sumSigns = 0;
    }

    BBarIntegrationPoint BBarIntegrator::generateIntegrationPoint(
        const BornPhSpPt& born, 
        const BornChannel& bornChannel, 
        const double amp2Born, 
        const double muF2,
        const double muR2,
        const std::array<double, 3>& unitCube
    ) const
    {
        BBarIntegrationPoint point = BBarIntegrationPoint{ born, bornChannel };

        point.u1 = unitCube[0];
        point.u2 = unitCube[1];
        point.u3 = unitCube[2];

        point.y = -1.0 + 2.0 * point.u2;
        point.xi = m_realPhaseSpace->xiMax(born, point.y) * point.u1;
        point.phi = 2.0 * PI * point.u3;

        point.unitCubeJacobian = 4.0 * PI * point.xi / point.u1;        

        point.zLeg1 = born.x1Bar + (1.0 - born.x1Bar) * point.u1;
        point.zLeg2 = born.x2Bar + (1.0 - born.x2Bar) * point.u1;

        point.f1Born = m_config.PDF->xfxQ2(bornChannel.id1, born.x1Bar, muF2) / born.x1Bar; 
        point.f2Born = m_config.PDF->xfxQ2(bornChannel.id2, born.x2Bar, muF2) / born.x2Bar; 

        // TODO probably these real luminosities are only for collinear remnants, not for the real part
        point.f1RealQ = m_config.PDF->xfxQ2(bornChannel.id1, born.x1Bar / point.zLeg1, muF2) / (born.x1Bar / point.zLeg1); 
        point.f2RealQ = m_config.PDF->xfxQ2(bornChannel.id2, born.x2Bar / point.zLeg2, muF2) / (born.x2Bar / point.zLeg2); 

        point.f1RealG = m_config.PDF->xfxQ2(21, born.x1Bar / point.zLeg1, muF2) / (born.x1Bar / point.zLeg1); 
        point.f2RealG = m_config.PDF->xfxQ2(21, born.x2Bar / point.zLeg2, muF2) / (born.x2Bar / point.zLeg2); 
        
        point.alphaS = m_config.alphaS(muR2);

        point.amp2Born = amp2Born;

        return point;
    }

} // namespace powheg_dy
