#include "bbar_integrator.h"

#include "powheg_dy/process.h"
#include "powheg_dy/math/rand.h"
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

    // TODO: Generalize and move
    std::vector<FKSRegion> drellYanRegions(const BornChannel& born, const RealChannel& real)  
    {
        std::vector<FKSRegion> regions;

        const int emittedIndex = static_cast<int>(real.outIDs.size()) - 1;
        const int emittedID = real.outIDs[emittedIndex];

        // q qbar -> X g
        if (real.id1 == born.id1
            && real.id2 == born.id2
            && real.outIDs.size() == born.outIDs.size() + 1
            && emittedID == 21)
        {
            regions.push_back({
                .fksPartonIdx = emittedIndex,   // if this indexes real.outIDs
                .emitter = 0,                // both ISR legs
                .hasSoft = true,
                .hasCollinearLeg1 = true,
                .hasCollinearLeg2 = true,
                .hasCollinearFinal = false
            });

            return regions;
        }

        // q g -> X q
        // leg 2: g -> born.id2 + emitted
        if (real.id1 == born.id1
            && real.id2 == 21
            && real.outIDs.size() == born.outIDs.size() + 1
            && emittedID == -born.id2)
        {
            regions.push_back({
                .fksPartonIdx = emittedIndex,
                .emitter = 2,
                .hasSoft = false,
                .hasCollinearLeg1 = false,
                .hasCollinearLeg2 = true,
                .hasCollinearFinal = false
            });

            return regions;
        }

        // g qbar -> X qbar
        // leg 1: g -> born.id1 + emitted
        if (real.id1 == 21
            && real.id2 == born.id2
            && real.outIDs.size() == born.outIDs.size() + 1
            && emittedID == -born.id1)
        {
            regions.push_back({
                .fksPartonIdx = emittedIndex,
                .emitter = 1,
                .hasSoft = false,
                .hasCollinearLeg1 = true,
                .hasCollinearLeg2 = false,
                .hasCollinearFinal = false
            });

            return regions;
        }

        throw std::runtime_error("Real channel does not match expected DY real channels.");
    }

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
        
        dSigma += m_bornVirtual.dSigmaVirtual(point, muF2, muR2);

        for (const auto& realChannel : m_process.realChannels(bornChannel))
        {
            for (const auto& fksRegion : drellYanRegions(bornChannel, realChannel))
            {
                dSigma += m_real.dSigmaRealMinusCT(point, realChannel, fksRegion, muF2, muR2, *m_realPhaseSpace);

                for (const auto& collinearChannel : remnantChannelsFromRegion(bornChannel, realChannel, fksRegion))
                    dSigma += m_collRemn.dSigmaCollinearRemnants(point, collinearChannel, muF2);
            }
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
        BBarIntegrationPoint point;

        point.u1 = unitCube[0];
        point.u2 = unitCube[1];
        point.u3 = unitCube[2];

        point.y = -1.0 + 2.0 * point.u2;
        point.xi = m_realPhaseSpace->xiMax(born, point.y) * point.u1;
        point.phi = 2.0 * PI * point.u3;

        point.born = born;
        point.bornChannel = bornChannel;
        point.real = m_realPhaseSpace->reconstruct(born, { point.xi, point.y, point.phi });

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
