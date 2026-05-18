#include "bbar_integrator.h"

#include "powheg/process.h"
#include "powheg/math/rand.h"
#include "powheg/bbar/collinear_remnant_channel.h"

namespace powheg
{
namespace 
{
    static constexpr double GEV_SQ_TO_PB = 0.389379338e9;
    static constexpr double SECURITY_FACTOR = 2.5;
    static constexpr int MAX_TRIALS = 10000;
    static constexpr double TINY_XI = 1.0e-6;
    static constexpr double TINY_Y = 1.0e-6;

} // anonymous namespace 
    
    BornEvent BBarIntegrator::computeWeightAndSampleChannel(const BornPhSpPt& born) const
    {
        const std::array<double, 3> unitCube = { rand(), rand(), rand() };

        std::vector<WeightedChannel> weights = {};
        const double totalAbsoluteWeight = bTilde(born, unitCube, weights);

        // Sample the parton channel by the relative contribution to dSigma
        double u = rand(0.0, totalAbsoluteWeight);
        for (const WeightedChannel& weightedChannel : weights)
        {
            if (u < weightedChannel.absoluteWeight)
                return { born, weightedChannel.channel, totalAbsoluteWeight, weightedChannel.weightSign };

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
                Log::warn << "Born weight " << bornEvent.absoluteWeight << " exceeds max weight " << m_maxWeight << ", accepting event." << Log::endl;

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
        const std::array<double, 3>& unitCube,
        std::vector<WeightedChannel>& weights
    ) const
    {
        const double muF2 = born.sHat;
        const double muR2 = born.sHat;

        // const auto& bornChannels = m_process.getBornChannelsOld();
        // std::vector<BBarCacheOld> caches;

        // caches.reserve(bornChannels.size());
        
        // prepare caches
        
        const BBarCache cache = prepareCache(unitCube, born);

        weights.reserve(cache.channels.size());

        double totalAbsoluteWeight = 0.0;
        for (size_t channelIdx = 0; channelIdx < cache.channels.size(); channelIdx++)
        {
            double channelWeight = 0.0;

            if (m_config.BTILDE_USE_BORN) channelWeight += m_bornVirtual.dSigmaBorn(cache, channelIdx);
            if (m_config.BTILDE_USE_VIRT) channelWeight += m_bornVirtual.dSigmaVirtual(cache, channelIdx);
            if (m_config.BTILDE_USE_REAL) channelWeight += m_real.dSigmaRealMinusCT(cache, channelIdx);
            if (m_config.BTILDE_USE_COLL) channelWeight += m_collRemn.dSigmaCollinearRemnants(cache, channelIdx);

            totalAbsoluteWeight += std::abs(channelWeight);
            weights.push_back({ cache.channels[channelIdx].channel.bornChannel, 
                cache.channels[channelIdx].amp2Born, channelWeight >= 0, std::abs(channelWeight) });
        }

        
        // for (size_t bornChannelIdx = 0; bornChannelIdx < bornChannels.size(); bornChannelIdx++)
        // {
        //     const BornChannel& bornChannel = bornChannels[bornChannelIdx];

        //     const double amp2Born = m_process.bornAmp2(born, bornChannel);
        //     const BBarCacheOld point = prepareCache(born, bornChannel, amp2Born, muF2, muR2, unitCube);

        //     double dSigma = m_bornVirtual.dSigmaBorn(point);

        //     if (!m_config.BTILDE_BORNONLY)
        //     {
        //         dSigma += m_bornVirtual.dSigmaVirtual(point, muF2, muR2);
        //         // dSigma += m_real.dSigmaRealMinusCT(point, bornChannelIdx, muF2);

        //         // const auto& realChannels = m_process.getRealChannelsOld();
        //         // for (size_t realChannelIdx = 0; realChannelIdx < realChannels.size(); realChannelIdx++)
        //         // {
        //             // const RealChannel& realChannel = realChannels[realChannelIdx];
        //             // for (const auto& fksRegion : m_process.getFKSRegionsOld(bornChannelIdx, realChannelIdx))
        //                 // dSigma += m_collRemn.dSigmaCollinearRemnants(point, fksRegion, muF2);
        //         // }
        //     }

        //     const double absoluteWeight = std::abs(dSigma);
        //     const int sign = dSigma >= 0;

        //     totalAbsoluteWeight += absoluteWeight;
        //     weights.push_back({ bornChannel, amp2Born, sign, absoluteWeight });
        // }

        return totalAbsoluteWeight;
    }

    RadiationCache BBarIntegrator::sampleRadiationVariables(
        const std::array<double, 3>& unitCube,
        const BornPhSpPt& born
    ) const
    {
        RadiationCache rad;

        double unitCubeJacobianTilde = 1.0;

        rad.xiTilde = (3.0 - 2.0 * unitCube[0]) * unitCube[0] * unitCube[0];
        rad.xiTilde = rad.xiTilde * (1.0 - 2.0 * TINY_XI) + TINY_XI;
        unitCubeJacobianTilde *= 6.0 * (1.0 - unitCube[0]) * unitCube[0];

        rad.xiMaxLeg1 = 1.0 - born.x1Bar;
        rad.xiMaxLeg2 = 1.0 - born.x2Bar;

        rad.xiLeg1 = rad.xiTilde * rad.xiMaxLeg1;
        rad.xiLeg2 = rad.xiTilde * rad.xiMaxLeg2;

        rad.zLeg1 = 1.0 - rad.xiLeg1;
        rad.zLeg2 = 1.0 - rad.xiLeg2;

        rad.dzduLeg1 = (6.0 * (1.0 - unitCube[0]) * unitCube[0]) * rad.xiMaxLeg1;
        rad.dzduLeg2 = (6.0 * (1.0 - unitCube[0]) * unitCube[0]) * rad.xiMaxLeg2;

        if (!m_config.BTILDE_USE_REAL)
            return rad;

        rad.y = 1.0 - 2.0 * unitCube[1];
        unitCubeJacobianTilde *= 2.0;
        unitCubeJacobianTilde *= 1.5 * (1.0 - rad.y * rad.y);
        rad.y = 1.5 * (rad.y - std::pow(rad.y, 3) / 3) * (1.0 - TINY_Y);

        rad.xiMax = m_realPhaseSpace.xiMax(born, rad.y);
        rad.xi = rad.xiTilde * rad.xiMax;

        rad.phi = 2.0 * PI * unitCube[2];
        unitCubeJacobianTilde *= 2.0 * PI;

        rad.real = m_realPhaseSpace.reconstruct(born, { rad.xi, rad.y, rad.phi });

        rad.jacobianRad = unitCubeJacobianTilde * rad.xiMax * rad.real.radJacobian;
        rad.jacobianRadSoftOverXiSq = unitCubeJacobianTilde * born.sHat / std::pow(4.0 * PI, 3);

        const double j1 = born.sHat / std::pow(4.0 * PI, 3) * rad.xiLeg1 / (1.0 - rad.xiLeg1);
        const double j2 = born.sHat / std::pow(4.0 * PI, 3) * rad.xiLeg2 / (1.0 - rad.xiLeg2);

        rad.jacobianRadLeg1 = unitCubeJacobianTilde * rad.xiMaxLeg1 * j1;
        rad.jacobianRadLeg2 = unitCubeJacobianTilde * rad.xiMaxLeg2 * j2;

        return rad;
    }

    ChannelCache BBarIntegrator::prepareChannelCache(
        const BBarCache& cache, 
        const Channel& channel
    ) const
    {
        ChannelCache channelCache{ channel };

        channelCache.amp2Born = m_process.bornAmp2(cache.born, channel.bornChannel);

        if (m_config.BTILDE_USE_VIRT)
            channelCache.amp2Virt = m_process.virtAmp2(cache.born, channel.bornChannel, channelCache.amp2Born, cache.muR2);

        channelCache.f1Born = m_config.PDF->xfxQ2(channel.bornChannel.id1, cache.born.x1Bar, cache.muF2) / cache.born.x1Bar; 
        channelCache.f2Born = m_config.PDF->xfxQ2(channel.bornChannel.id2, cache.born.x2Bar, cache.muF2) / cache.born.x2Bar; 

        channelCache.luminosityBorn = channelCache.f1Born * channelCache.f2Born;

        if (m_config.BTILDE_USE_REAL || m_config.BTILDE_USE_COLL)
        {
            channelCache.f1Coll = std::vector<double>(channel.relevantRealChannels.size());
            channelCache.f2Coll = std::vector<double>(channel.relevantRealChannels.size());

            for (size_t realIdx = 0; realIdx < channel.relevantRealChannels.size(); realIdx++)
            {
                const RealChannel& realChannel = channel.relevantRealChannels[realIdx];

                if (channel.hasCollinearLeg1[realIdx])
                {
                    const double x1Leg1 = cache.born.x1Bar / cache.rad.zLeg1;
                    channelCache.f1Coll[realIdx] = m_config.PDF->xfxQ2(realChannel.id1, x1Leg1, cache.muF2) / x1Leg1;
                }
                if (channel.hasCollinearLeg2[realIdx])
                {
                    const double x2Leg2 = cache.born.x2Bar / cache.rad.zLeg2;
                    channelCache.f2Coll[realIdx] = m_config.PDF->xfxQ2(realChannel.id2, x2Leg2, cache.muF2) / x2Leg2;
                }
            }
        }

        return channelCache;
    }

    BBarCache BBarIntegrator::prepareCache(
        const std::array<double, 3>& unitCube, 
        const BornPhSpPt& born
    ) const
    {
        BBarCache cache{ born };

        cache.muF2 = born.sHat;
        cache.muR2 = born.sHat;
        cache.alphaS = m_config.alphaS(cache.muR2);

        cache.fluxFactorBorn = 1.0 / (2.0 * born.sHat);

        if (m_config.BTILDE_USE_REAL || m_config.BTILDE_USE_COLL)
            cache.rad = sampleRadiationVariables(unitCube, born);

        const std::vector<Channel> channels = m_process.getChannels();
        cache.channels.reserve(channels.size());

        for (const Channel& channel : channels)
            cache.channels.push_back(prepareChannelCache(cache, channel));

        return cache;
    }

    void BBarIntegrator::determineMaxWeight()
    {
        if (m_config.BORN_VETO_WEIGHT != -1.0)
        {
            m_maxWeight = m_config.BORN_VETO_WEIGHT;
            Log::info << "Using manual Born Veto weight " << m_maxWeight << Log::endl;
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

            Log::info << "Done determining Born Veto weight." << Log::endl << Log::endl;
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

    // BBarCacheOld BBarIntegrator::prepareCache(
    //     const BornPhSpPt& born, 
    //     const BornChannel& bornChannel, 
    //     const double amp2Born, 
    //     const double muF2,
    //     const double muR2,
    //     const std::array<double, 3>& unitCube
    // ) const
    // {
    //     BBarCacheOld cache;

    //     cache.born = born;
    //     cache.bornChannel = bornChannel;

    //     cache.alphaS = m_config.alphaS(muR2);
    //     cache.amp2Born = amp2Born;

    //     const double TINY_XI = 1.0e-6;
    //     const double TINY_Y = 1.0e-6;
        
    //     ////////////////////////////////// For the Counterterms //////////////////////////////////

    //     double unitCubeJacobianTilde = 1.0;

    //     cache.xiTilde = (3.0 - 2.0 * unitCube[0]) * unitCube[0] * unitCube[0];
    //     cache.xiTilde = cache.xiTilde * (1.0 - 2.0 * TINY_XI) + TINY_XI;
    //     unitCubeJacobianTilde *= 6.0 * (1.0 - unitCube[0]) * unitCube[0];

    //     cache.y = 1.0 - 2.0 * unitCube[1];
    //     unitCubeJacobianTilde *= 2.0;
    //     unitCubeJacobianTilde *= 1.5 * (1.0 - cache.y * cache.y);
    //     cache.y = 1.5 * (cache.y - std::pow(cache.y, 3) / 3) * (1.0 - TINY_Y);

    //     cache.xiMax = m_realPhaseSpace.xiMax(cache.born, cache.y);
    //     cache.xiMaxLeg1 = 1.0 - cache.born.x1Bar;
    //     cache.xiMaxLeg2 = 1.0 - cache.born.x2Bar;

    //     cache.xi = cache.xiTilde * cache.xiMax;
    //     cache.xiLeg1 = cache.xiTilde * cache.xiMaxLeg1;
    //     cache.xiLeg2 = cache.xiTilde * cache.xiMaxLeg2;

    //     cache.zLeg1 = 1.0 - cache.xiLeg1;
    //     cache.zLeg2 = 1.0 - cache.xiLeg2;

    //     cache.dzduLeg1 = (6.0 * (1.0 - unitCube[0]) * unitCube[0]) * cache.xiMaxLeg1;
    //     cache.dzduLeg2 = (6.0 * (1.0 - unitCube[0]) * unitCube[0]) * cache.xiMaxLeg2;

    //     cache.phi = 2.0 * PI * unitCube[2];
    //     unitCubeJacobianTilde *= 2.0 * PI;

    //     cache.real = m_realPhaseSpace.reconstruct(born, { cache.xi, cache.y, cache.phi });

    //     cache.jacobian = cache.born.jacobian * unitCubeJacobianTilde * cache.xiMax * cache.real.radJacobian;
    //     cache.jacobianOverXiSqSoft = cache.born.jacobian * unitCubeJacobianTilde * cache.born.sHat / std::pow(4.0 * PI, 3);

    //     const double jacobianRadLeg1 = cache.born.sHat / std::pow(4.0 * PI, 3) * cache.xiLeg1 / (1.0 - cache.xiLeg1);
    //     const double jacobianRadLeg2 = cache.born.sHat / std::pow(4.0 * PI, 3) * cache.xiLeg2 / (1.0 - cache.xiLeg2);

    //     cache.jacobianLeg1 = unitCubeJacobianTilde * cache.xiMaxLeg1 * cache.born.jacobian * jacobianRadLeg1;
    //     cache.jacobianLeg2 = unitCubeJacobianTilde * cache.xiMaxLeg2 * cache.born.jacobian * jacobianRadLeg2;

    //     /////////////////////////////// For the Collinear Remnants ///////////////////////////////

    //     cache.f1Born = m_config.PDF->xfxQ2(bornChannel.id1, born.x1Bar, muF2) / born.x1Bar; 
    //     cache.f2Born = m_config.PDF->xfxQ2(bornChannel.id2, born.x2Bar, muF2) / born.x2Bar; 

    //     cache.f1RealQ = m_config.PDF->xfxQ2(bornChannel.id1, born.x1Bar / cache.zLeg1, muF2) / (born.x1Bar / cache.zLeg1); 
    //     cache.f2RealQ = m_config.PDF->xfxQ2(bornChannel.id2, born.x2Bar / cache.zLeg2, muF2) / (born.x2Bar / cache.zLeg2); 

    //     cache.f1RealG = m_config.PDF->xfxQ2(21, born.x1Bar / cache.zLeg1, muF2) / (born.x1Bar / cache.zLeg1); 
    //     cache.f2RealG = m_config.PDF->xfxQ2(21, born.x2Bar / cache.zLeg2, muF2) / (born.x2Bar / cache.zLeg2); 

    //     return cache;
    // }

} // namespace powheg
