#include "real_minus_ct_handler.h"

#include "powheg/process.h"

namespace powheg
{
    double RealMinusCTHandler::dSigmaRealMinusCT(
        const BBarCache& cache,
        size_t channelIdx
    ) const
    {
        const ChannelCache& channelCache = cache.channels[channelIdx];

        double dSigma = 0.0;
        
        const auto& realChannels = channelCache.channel.relevantRealChannels;
        for (size_t realChannelIdx = 0; realChannelIdx < realChannels.size(); realChannelIdx++)
        {
            const RealChannel& realChannel = realChannels[realChannelIdx];

            const auto& regions = channelCache.channel.fksRegions[realChannelIdx];

            // TODO: If a real channel is associated with multiple born channels, this gets computed multiple times
            const double realAmp2 = m_process.realAmp2(cache.rad.real, realChannels[realChannelIdx], cache.alphaS);
            const double f1 = m_config.PDF->xfxQ2(realChannel.id1, cache.rad.real.x1, cache.muF2) / cache.rad.real.x1;
            const double f2 = m_config.PDF->xfxQ2(realChannel.id2, cache.rad.real.x2, cache.muF2) / cache.rad.real.x2;
            const double fluxFactor = 1.0 / (2.0 * cache.rad.real.sHatReal);

            for (const FKSRegion& region : regions)
            {
                if (region.fsrRegion.has_value())
                    throw std::runtime_error("FSR not yet implemented.");

                dSigma += fksPartition(region, cache.rad.y) * cache.rad.jacobianRad * cache.born.jacobian * f1 * f2 * fluxFactor * realAmp2;

                if (region.soft)
                {
                    const double softCT = softCounterterm(cache, channelIdx, region);
                    const double boundaryTerm = softCT * std::log(cache.rad.xiMax) * cache.rad.xiTilde;

                    dSigma += -softCT + boundaryTerm;
                }

                if (region.isrRegion.has_value() && region.isrRegion->leg == Leg::LEG1)
                {
                    const double leg1CollinearCT = leg1CollinearCounterterm(cache, channelIdx, realChannelIdx, region);
                    dSigma += -leg1CollinearCT;

                    if (region.soft)
                    {
                        const double leg1SoftCollinearCT = leg1SoftCollinearCounterterm(cache, channelIdx, region);
                        const double boundaryTerm = -leg1SoftCollinearCT * std::log(cache.rad.xiMaxLeg1) * cache.rad.xiTilde;
                        dSigma += leg1SoftCollinearCT + boundaryTerm;
                    }
                }

                if (region.isrRegion.has_value() && region.isrRegion->leg == Leg::LEG2)
                {
                    const double leg2CollinearCT = leg2CollinearCounterterm(cache, channelIdx, realChannelIdx, region);
                    dSigma += -leg2CollinearCT;

                    if (region.soft)
                    {
                        const double leg2SoftCollinearCT = leg2SoftCollinearCounterterm(cache, channelIdx, region);
                        const double boundaryTerm = -leg2SoftCollinearCT * std::log(cache.rad.xiMaxLeg2) * cache.rad.xiTilde;
                        dSigma += leg2SoftCollinearCT + boundaryTerm;
                    }
                }
            }
        }

        return dSigma;
    }

    double RealMinusCTHandler::dSigmaReal(
        const BBarCacheOld& cache,
        const RealChannel& realChannel,
        const double muF2
    ) const
    {
        const double amp2 = m_process.realAmp2(cache.real, realChannel, cache.alphaS);
        const double fluxFactor = 1 / (2.0 * cache.real.sHatReal);

        const double f1 = m_config.PDF->xfxQ2(realChannel.id1, cache.real.x1, muF2) / cache.real.x1;
        const double f2 = m_config.PDF->xfxQ2(realChannel.id2, cache.real.x2, muF2) / cache.real.x2;
        const double luminosity = f1 * f2;

        return cache.jacobian * luminosity * fluxFactor * amp2;
    }

    double RealMinusCTHandler::softCounterterm(
        const BBarCache& cache,
        size_t channelIdx,
        const FKSRegion& region
    ) const
    {
        const ChannelCache& channelCache = cache.channels[channelIdx];

        const FourVector softGluonVec = {
            1.0,
            std::sqrt(1 - cache.rad.y * cache.rad.y) * std::sin(cache.rad.phi),
            std::sqrt(1 - cache.rad.y * cache.rad.y) * std::cos(cache.rad.phi),
            cache.rad.y
        };

        // TODO: Generalize from here
        
        const double B12 = m_config.C_F * cache.fluxFactorBorn * channelCache.amp2Born;

        const double eikonalFactor = 2.0 * B12 
            * dot(cache.born.p1Bar, cache.born.p2Bar)
            / dot(cache.born.p1Bar, softGluonVec)
            / dot(cache.born.p2Bar, softGluonVec);

        const double fluxFactor = 1.0 / (2.0 * cache.born.sHat);
        const double luminosity = channelCache.f1Born * channelCache.f2Born;

        const double coupling = 4.0 * PI * cache.alphaS;
        const double amp2SoftTimesXiSq = 8.0 * coupling * eikonalFactor;

        // this returns rrrs in sigreal.f
        return fksPartition(region, cache.rad.y) * cache.born.jacobian * cache.rad.jacobianRadSoftOverXiSq 
            * luminosity * fluxFactor * amp2SoftTimesXiSq;
    }

    double RealMinusCTHandler::leg1CollinearCounterterm(
        const BBarCache& cache,
        size_t channelIdx,
        size_t realChannelIdx,
        const FKSRegion& region
    ) const
    {
        const ChannelCache& channelCache = cache.channels[channelIdx];

        const double z  = cache.rad.zLeg1;
        const double x1 = cache.born.x1Bar / z;

        // const double f1 = m_config.PDF->xfxQ2(realChannel.id1, x1, muF2) / x1;
        const double luminosity = channelCache.f1Coll[realChannelIdx] * channelCache.f2Born;

        double kernel = 0.0;
        if (region.isrRegion->splitting == ISRSplitting::QQ)
            kernel = oneMinusZTimesPqq(z);
        else if (region.isrRegion->splitting == ISRSplitting::GQ)
            kernel = oneMinusZTimesPgq(z);
        else
            throw std::runtime_error("This channel has no leg 1 collinear singularity");

        const double fluxFactor = 1.0 / (2.0 * cache.born.sHat);

        const double coupling = 4.0 * PI * cache.alphaS;
        const double singularFactor = 1.0 / cache.rad.xiLeg1 / cache.rad.xiLeg1 / (1.0 - cache.rad.y);

        const double amp2Collinear =
            4.0 / cache.born.sHat
            * coupling
            * kernel
            * singularFactor
            * channelCache.amp2Born;

        // this returns rrrp on line 132 of sigreal.f
        return cache.born.jacobian * cache.rad.jacobianRadLeg1 * luminosity * fluxFactor * amp2Collinear;
    }

    double RealMinusCTHandler::leg2CollinearCounterterm(
        const BBarCache& cache,
        size_t channelIdx,
        size_t realChannelIdx,
        const FKSRegion& region
    ) const
    {
        const ChannelCache& channelCache = cache.channels[channelIdx];

        const double z  = cache.rad.zLeg2;
        const double x2 = cache.born.x2Bar / z;

        const double luminosity = channelCache.f1Born * channelCache.f2Coll[realChannelIdx];

        // const FourVector kPerp = { 0.0, std::sin(cache.phi), std::cos(cache.phi), 0.0 };

        double kernel = 0.0;
        if (region.isrRegion->splitting == ISRSplitting::QQ)
            kernel = oneMinusZTimesPqq(z);
        else if (region.isrRegion->splitting == ISRSplitting::GQ)
            kernel = oneMinusZTimesPgq(z);
        else
            throw std::runtime_error("This channel has no leg 2 collinear singularity");

        const double fluxFactor = 1.0 / (2.0 * cache.born.sHat);

        const double coupling = 4.0 * PI * cache.alphaS;
        const double singularFactor = 1.0 / cache.rad.xiLeg2 / cache.rad.xiLeg2 / (1.0 + cache.rad.y);

        const double amp2Collinear =
            4.0 / cache.born.sHat
            * coupling
            * kernel
            * singularFactor
            * channelCache.amp2Born;

        // this returns rrrm on line 142 of sigreal.f
        return cache.born.jacobian * cache.rad.jacobianRadLeg2 * luminosity * fluxFactor * amp2Collinear;
    }

    double RealMinusCTHandler::leg1SoftCollinearCounterterm(
        const BBarCache& cache,
        size_t channelIdx,
        const FKSRegion& region
    ) const
    {
        const ChannelCache& channelCache = cache.channels[channelIdx];

        double kernel = 0.0;
        if (region.isrRegion->splitting == ISRSplitting::QQ)
            kernel = oneMinusZTimesPqqAtOne();
        else if (region.isrRegion->splitting == ISRSplitting::GQ)
            kernel = oneMinusZTimesPgqAtOne();
        else
            throw std::runtime_error("This channel has no leg 1 soft collinear singularity");

        const double luminosity = channelCache.f1Born * channelCache.f2Born;

        const double coupling = 4.0 * PI * cache.alphaS;
        const double singularFactorTimesXiSq = 1.0 / (1.0 - cache.rad.y);

        const double amp2SoftCollinearTimesXiSq =
            4.0 / cache.born.sHat
            * coupling
            * kernel
            * singularFactorTimesXiSq
            * channelCache.amp2Born;
            
        // this returns rrrps on line 134 of sigreal.f
        return cache.born.jacobian * cache.rad.jacobianRadSoftOverXiSq * luminosity * cache.fluxFactorBorn * amp2SoftCollinearTimesXiSq;
    }

    double RealMinusCTHandler::leg2SoftCollinearCounterterm(
        const BBarCache& cache,
        size_t channelIdx,
        const FKSRegion& region
    ) const
    {
        const ChannelCache& channelCache = cache.channels[channelIdx];

        double kernel = 0.0;
        if (region.isrRegion->splitting == ISRSplitting::QQ)
            kernel = oneMinusZTimesPqqAtOne();
        else if (region.isrRegion->splitting == ISRSplitting::GQ)
            kernel = oneMinusZTimesPgqAtOne();
        else
            throw std::runtime_error("This channel has no leg 2 soft collinear singularity");

        const double fluxFactor = 1.0 / (2.0 * cache.born.sHat);
        const double luminosity = channelCache.f1Born * channelCache.f2Born;

        const double coupling = 4.0 * PI * cache.alphaS;
        const double singularFactorTimesXiSq = 1.0 / (1.0 + cache.rad.y);

        const double amp2SoftCollinearTimesXiSq =
            4.0 / cache.born.sHat
            * coupling
            * kernel
            * singularFactorTimesXiSq
            * channelCache.amp2Born;
            
        // this returns rrrms on line 144 of sigreal.f
        return cache.born.jacobian * cache.rad.jacobianRadSoftOverXiSq * luminosity * fluxFactor * amp2SoftCollinearTimesXiSq;
    }
    
    double RealMinusCTHandler::fksPartition(
        const FKSRegion& region,
        double y
    ) const
    {
        if (!region.soft)
            return 1.0;

        if (region.fsrRegion.has_value())
            throw std::runtime_error("FSR not yet implemented");

        if (region.isrRegion->leg == Leg::LEG1)
            return 0.5 * (1.0 + y);

        if (region.isrRegion->leg == Leg::LEG2)
            return 0.5 * (1.0 - y);

        throw std::runtime_error("Invalid ISR leg.");
    }

    double RealMinusCTHandler::oneMinusZTimesPqq(const double z) const
    {
        // eq. (2.103) times (1 - z)
        return m_config.C_F * (1.0 + z * z);
    }

    double RealMinusCTHandler::oneMinusZTimesPqg(const double z) const
    {
        // eq. (2.104) times (1 - z)
        return m_config.C_F * (1.0 + (1.0 - z) * (1.0 - z)) * (1.0 - z) / z;
    }

    double RealMinusCTHandler::oneMinusZTimesPgq(const double z) const
    {
        // eq. (2.105) times (1 - z)
        return m_config.T_F * (1.0 - 2.0 * z * (1.0 - z)) * (1.0 - z);
    }

    double RealMinusCTHandler::oneMinusZTimesPgg(const double z) const
    {
        // eq. (2.106) times (1 - z)
        return 2.0 * m_config.C_A * (z + (1.0 - z) * (1.0 - z) / z + z * (1.0 - z) * (1.0 - z));
    }

    double RealMinusCTHandler::oneMinusZTimesPqqAtOne() const
    {
        // eq. (2.103) times (1 - z) at z = 1
        return 2.0 * m_config.C_F;;
    }

    double RealMinusCTHandler::oneMinusZTimesPqgAtOne() const
    {
        // eq. (2.104) times (1 - z) at z = 1
        return 0.0;
    }

    double RealMinusCTHandler::oneMinusZTimesPgqAtOne() const
    {
        // eq. (2.105) times (1 - z) at z = 1
        return 0.0;
    }

    double RealMinusCTHandler::oneMinusZTimesPggAtOne() const
    {
        // eq. (2.106) times (1 - z) at z = 1
        return 2.0 * m_config.C_A;
    }
   
} // namespace powheg
