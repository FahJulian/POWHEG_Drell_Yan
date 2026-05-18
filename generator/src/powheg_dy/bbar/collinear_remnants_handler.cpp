#include "collinear_remnants_handler.h"

#include "powheg_dy/process.h"

namespace powheg
{
    double CollinearRemnantsHandler::dSigmaCollinearRemnants(
        const BBarIntegrationPoint& point,
        const CollinearRemnantChannel& channel,
        const double muF2
    ) const
    {
        if (!m_config.BTILDE_USE_COLL_REMNANTS)
            return 0.0;

        if (channel.leg == 1)
            return dSigmaCollinearRemnantsLeg1(point, channel.splitting, muF2);
        else if (channel.leg == 2)
            return dSigmaCollinearRemnantsLeg2(point, channel.splitting, muF2);
            
        throw std::runtime_error("Invalid collinear remnant leg");
        return 0.0;
    }

    double CollinearRemnantsHandler::dSigmaCollinearRemnantsLeg1(
        const BBarIntegrationPoint& point,
        const CollinearRemnantSplitting& splitting,
        const double muF2
    ) const
    {
        const double z = point.zLeg1CollRemn;
        const double s = point.born.sHat / z;
        const double dzdu = 1.0 - point.born.x1Bar;

        const double dSigmaBornPartonic = 1.0 / (2.0 * point.born.sHat) * point.born.jacobian * point.amp2Born;
        const double luminosityOne = point.f1Born * point.f2Born;
        
        double luminosityZ;
        double oneMinusZTimesPOne;
        double oneMinusZTimesPZ;
        double derivativePZ;

        if (splitting == CollinearRemnantSplitting::QQ)
        {
            luminosityZ        = point.f1RealQ * point.f2Born;
            oneMinusZTimesPOne = oneMinusZTimesPqq(1.0);
            oneMinusZTimesPZ   = oneMinusZTimesPqq(z);
            derivativePZ       = derivativePqq(z);
        }
        else if (splitting == CollinearRemnantSplitting::QG)
        {
            luminosityZ        = point.f1RealQ * point.f2Born;
            oneMinusZTimesPOne = oneMinusZTimesPqg(1.0);
            oneMinusZTimesPZ   = oneMinusZTimesPqg(z);
            derivativePZ       = derivativePqg(z);
        }
        else if (splitting == CollinearRemnantSplitting::GQ)
        {
            luminosityZ        = point.f1RealG * point.f2Born;
            oneMinusZTimesPOne = oneMinusZTimesPgq(1.0);
            oneMinusZTimesPZ   = oneMinusZTimesPgq(z);
            derivativePZ       = derivativePgq(z);
        }
        else if (splitting == CollinearRemnantSplitting::GG)
        {
            luminosityZ        = point.f1RealG * point.f2Born;
            oneMinusZTimesPOne = oneMinusZTimesPgg(1.0);
            oneMinusZTimesPZ   = oneMinusZTimesPgg(z);
            derivativePZ       = derivativePgg(z);
        }
        else
            throw std::runtime_error("Invalid Collinear Remnant Splitting");

        const double logS = std::log(s / muF2);
        const double logSHat = std::log(point.born.sHat / muF2);

        const double fZ = luminosityZ * oneMinusZTimesPZ / z
            * (logS / (1.0 - z) + 2.0 * std::log(1.0 - z) / (1.0 - z));
        const double f1 = luminosityOne * oneMinusZTimesPOne
            * (logSHat / (1.0 - z) + 2.0 * std::log(1.0 - z) / (1.0 - z));
        
        const double regular = - luminosityZ * derivativePZ / z;
        
        const double logZEndpoint = std::log(1.0 - point.born.x1Bar);
        const double endpoint = luminosityOne * oneMinusZTimesPOne 
            * (logSHat * logZEndpoint + logZEndpoint * logZEndpoint);

        return point.alphaS / (2.0 * PI) 
            * (dzdu * (fZ - f1 + regular) + endpoint)
            * dSigmaBornPartonic;
    }

    double CollinearRemnantsHandler::dSigmaCollinearRemnantsLeg2(
        const BBarIntegrationPoint& point,
        const CollinearRemnantSplitting& splitting,
        const double muF2
    ) const
    {
        const double z = point.zLeg2CollRemn;
        const double s = point.born.sHat / z;
        const double dzdu = 1.0 - point.born.x2Bar;

        const double dSigmaBornPartonic = 1.0 / (2.0 * point.born.sHat) * point.born.jacobian * point.amp2Born;
        const double luminosityOne = point.f1Born * point.f2Born;
        
        double luminosityZ;
        double oneMinusZTimesPOne;
        double oneMinusZTimesPZ;
        double derivativePZ;

        if (splitting == CollinearRemnantSplitting::QQ)
        {
            luminosityZ        = point.f1Born * point.f2RealQ;
            oneMinusZTimesPOne = oneMinusZTimesPqq(1.0);
            oneMinusZTimesPZ   = oneMinusZTimesPqq(z);
            derivativePZ       = derivativePqq(z);
        }
        else if (splitting == CollinearRemnantSplitting::QG)
        {
            luminosityZ        = point.f1Born * point.f2RealQ;
            oneMinusZTimesPOne = oneMinusZTimesPqg(1.0);
            oneMinusZTimesPZ   = oneMinusZTimesPqg(z);
            derivativePZ       = derivativePqg(z);
        }
        else if (splitting == CollinearRemnantSplitting::GQ)
        {
            luminosityZ        = point.f1Born * point.f2RealG;
            oneMinusZTimesPOne = oneMinusZTimesPgq(1.0);
            oneMinusZTimesPZ   = oneMinusZTimesPgq(z);
            derivativePZ       = derivativePgq(z);
        }
        else if (splitting == CollinearRemnantSplitting::GG)
        {
            luminosityZ        = point.f1Born * point.f2RealG;
            oneMinusZTimesPOne = oneMinusZTimesPgg(1.0);
            oneMinusZTimesPZ   = oneMinusZTimesPgg(z);
            derivativePZ       = derivativePgg(z);
        }
        else
            throw std::runtime_error("Invalid Collinear Remnant Splitting");

        const double logS = std::log(s / muF2);
        const double logSHat = std::log(point.born.sHat / muF2);

        const double fZ = luminosityZ * oneMinusZTimesPZ / z
            * (logS / (1.0 - z) + 2.0 * std::log(1.0 - z) / (1.0 - z));
        const double f1 = luminosityOne * oneMinusZTimesPOne
            * (logSHat / (1.0 - z) + 2.0 * std::log(1.0 - z) / (1.0 - z));
        
        const double regular = - luminosityZ * derivativePZ / z;
        
        const double logZEndpoint = std::log(1.0 - point.born.x2Bar);
        const double endpoint = luminosityOne * oneMinusZTimesPOne 
            * (logSHat * logZEndpoint + logZEndpoint * logZEndpoint);

        return point.alphaS / (2.0 * PI) 
            * (dzdu * (fZ - f1 + regular) + endpoint)
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

} // namespace powheg_dy
