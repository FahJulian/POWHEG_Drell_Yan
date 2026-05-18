#include "collinear_remnants_handler.h"

#include "powheg/process.h"

namespace powheg
{
    double CollinearRemnantsHandler::dSigmaCollinearRemnants(
        const BBarCache& cache,
        size_t channelIdx
    ) const
    {
        double dSigma = 0.0;
        
        const ChannelCache& channelCache = cache.channels[channelIdx];
        
        const auto& realChannels = channelCache.channel.relevantRealChannels;
        for (size_t realChannelIdx = 0; realChannelIdx < realChannels.size(); realChannelIdx++)
        {
            for (const FKSRegion& region : channelCache.channel.fksRegions[realChannelIdx])
            {
                if (!region.isrRegion.has_value())
                    continue;

                if (region.isrRegion->leg == Leg::LEG1)
                    dSigma += dSigmaCollinearRemnantsLeg1(region.isrRegion->splitting, cache, channelIdx, realChannelIdx);
                else if (region.isrRegion->leg == Leg::LEG2)
                    dSigma += dSigmaCollinearRemnantsLeg2(region.isrRegion->splitting, cache, channelIdx, realChannelIdx);
                else
                    throw std::runtime_error("Invalid collinear remnant leg");
            }
        }

        return dSigma;
    }

    double CollinearRemnantsHandler::dSigmaCollinearRemnantsLeg1(
        const ISRSplitting& splitting,
        const BBarCache& cache, 
        size_t channelIdx,
        size_t realChannelIdx
    ) const
    {
        const ChannelCache& channel = cache.channels[channelIdx];

        const double z = cache.rad.zLeg1;
        const double s = cache.born.sHat / z;

        const double dSigmaBornPartonic = cache.fluxFactorBorn * cache.born.jacobian * channel.amp2Born;
        const double luminosityOne = channel.f1Born * channel.f2Born;
        
        double luminosityZ;
        double oneMinusZTimesPOne;
        double oneMinusZTimesPZ;
        double derivativePZ;

        if (splitting == ISRSplitting::QQ)
        {
            luminosityZ        = channel.f1Coll[realChannelIdx] * channel.f2Born;
            oneMinusZTimesPZ   = m_config.C_F * (1.0 + z * z);
            oneMinusZTimesPOne = 2.0 * m_config.C_F;
            derivativePZ       = derivativePqq(z);
        }
        else if (splitting == ISRSplitting::GQ)
        {
            luminosityZ        = channel.f1Coll[realChannelIdx] * channel.f2Born;
            oneMinusZTimesPZ   = m_config.T_F * (1.0 - 2.0 * z * (1.0 - z)) * (1.0 - z);
            oneMinusZTimesPOne = 0.0;
            derivativePZ       = -m_config.T_F * 2.0 * z * (1.0 - z);
        }
        else if (splitting == ISRSplitting::QG)
            throw std::runtime_error("QG splittings not implemented");
        else if (splitting == ISRSplitting::GG)
            throw std::runtime_error("GG splittings not implemented");
        else
            throw std::runtime_error("Invalid Collinear Remnant Splitting");

        const double logS = std::log(s / cache.muF2);
        const double logSHat = std::log(cache.born.sHat / cache.muF2);

        // the distribution coefficient at z and z=1 in eq. (2.102)
        const double bracketZ = (logS + 2.0 * std::log(1.0 - z)) / (1.0 - z);
        const double bracketOne = (logSHat + 2.0 * std::log(1.0 - z)) / (1.0 - z);

        // the full distribution term in eq. (2.102)
        const double distribution = luminosityZ * oneMinusZTimesPZ * bracketZ / z
            - luminosityOne * oneMinusZTimesPOne * bracketOne;
        
        // the regular term in eq. (2.102)
        const double regular = - luminosityZ * derivativePZ / z;
        
        // because the z-integral is over z\in[x1bar, 1], there is a boundary term
        const double logZEndpoint = std::log(1.0 - cache.born.x1Bar);
        const double endpoint = luminosityOne * oneMinusZTimesPOne 
            * (logSHat * logZEndpoint + logZEndpoint * logZEndpoint);

        return cache.alphaS / (2.0 * PI)
            * (cache.rad.dzduLeg1 * (distribution + regular) + endpoint)
            * dSigmaBornPartonic;
    }

    double CollinearRemnantsHandler::dSigmaCollinearRemnantsLeg2(
        const ISRSplitting& splitting,
        const BBarCache& cache, 
        size_t channelIdx,
        size_t realChannelIdx
    ) const
    {
        const ChannelCache& channel = cache.channels[channelIdx];

        const double z = cache.rad.zLeg2;
        const double s = cache.born.sHat / z;

        const double dSigmaBornPartonic = cache.fluxFactorBorn * cache.born.jacobian * channel.amp2Born;
        const double luminosityOne = channel.f1Born * channel.f2Born;
        
        double luminosityZ;
        double oneMinusZTimesPOne;
        double oneMinusZTimesPZ;
        double derivativePZ;

        if (splitting == ISRSplitting::QQ)
        {
            luminosityZ        = channel.f1Born * channel.f2Coll[realChannelIdx];
            oneMinusZTimesPZ   = m_config.C_F * (1.0 + z * z);
            oneMinusZTimesPOne = 2.0 * m_config.C_F;
            derivativePZ       = derivativePqq(z);
        }
        else if (splitting == ISRSplitting::GQ)
        {
            luminosityZ        = channel.f1Born * channel.f2Coll[realChannelIdx];
            oneMinusZTimesPZ   = m_config.T_F * (1.0 - 2.0 * z * (1.0 - z)) * (1.0 - z);
            oneMinusZTimesPOne = 0.0;
            derivativePZ       = -m_config.T_F * 2.0 * z * (1.0 - z);
        }
        else if (splitting == ISRSplitting::QG)
            throw std::runtime_error("QG splittings not implemented");
        else if (splitting == ISRSplitting::GG)
            throw std::runtime_error("GG splittings not implemented");
        else
            throw std::runtime_error("Invalid Collinear Remnant Splitting");

        const double logS = std::log(s / cache.muF2);
        const double logSHat = std::log(cache.born.sHat / cache.muF2);

        // the distribution coefficient at z and z=1 in eq. (2.102)
        const double bracketZ = (logS + 2.0 * std::log(1.0 - z)) / (1.0 - z);
        const double bracketOne = (logSHat + 2.0 * std::log(1.0 - z)) / (1.0 - z);

        // the full distribution term in eq. (2.102)
        const double distribution = luminosityZ * oneMinusZTimesPZ * bracketZ / z
            - luminosityOne * oneMinusZTimesPOne * bracketOne;
        
        // the regular term in eq. (2.102)
        const double regular = - luminosityZ * derivativePZ / z;
        
        // because the z-integral is over z\in[x1bar, 1], there is a boundary term
        const double logZEndpoint = std::log(1.0 - cache.born.x2Bar);
        const double endpoint = luminosityOne * oneMinusZTimesPOne 
            * (logSHat * logZEndpoint + logZEndpoint * logZEndpoint);

        return cache.alphaS / (2.0 * PI)
            * (cache.rad.dzduLeg1 * (distribution + regular) + endpoint)
            * dSigmaBornPartonic;
    }

    double CollinearRemnantsHandler::oneMinusZTimesPqq(const double z) const
    {
        // eq. (2.103) times (1 - z)
        return m_config.C_F * (1.0 + z * z);
    }

    double CollinearRemnantsHandler::oneMinusZTimesPqg(const double z) const
    {
        // eq. (2.104) times (1 - z)
        return m_config.C_F * (1.0 + (1.0 - z) * (1.0 - z)) * (1.0 - z) / z;
    }

    double CollinearRemnantsHandler::oneMinusZTimesPgq(const double z) const
    {
        // eq. (2.105) times (1 - z)
        return m_config.T_F * (1.0 - 2.0 * z * (1.0 - z)) * (1.0 - z);
    }

    double CollinearRemnantsHandler::oneMinusZTimesPgg(const double z) const
    {
        // eq. (2.106) times (1 - z)
        return 2.0 * m_config.C_A * (z + (1.0 - z) * (1.0 - z) / z + z * (1.0 - z) * (1.0 - z));
    }

    double CollinearRemnantsHandler::derivativePqq(const double z) const
    {
        return - m_config.C_F * (1.0 - z);
    }

    double CollinearRemnantsHandler::derivativePqg(const double z) const
    {
        return - m_config.C_F * z;
    }

    double CollinearRemnantsHandler::derivativePgq(const double z) const
    {
        return - m_config.T_F * 2.0 * z * (1.0 - z);
    }

    double CollinearRemnantsHandler::derivativePgg(const double z) const
    {
        return 0.0;
    }

} // namespace powheg
