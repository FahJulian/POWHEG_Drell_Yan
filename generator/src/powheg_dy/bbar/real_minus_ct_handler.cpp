#include "real_minus_ct_handler.h"

#include "powheg_dy/process.h"

namespace powheg
{
namespace
{
    constexpr bool TEST = true;

} // anonymous namespace

    double RealMinusCTHandler::dSigmaRealMinusCT(
        const BBarIntegrationPoint& point,
        const RealChannel& realChannel,
        const FKSRegion& region,
        const double muF2
    ) const
    {
        if (region.hasCollinearFinal)
            throw std::runtime_error("FSR not yet implemented.");

        double dSigma = 0.0;
        
        if (m_config.BTILDE_USE_REAL)
            dSigma += dSigmaReal(point, realChannel, muF2);
        
        if (!m_config.BTILDE_USE_COUNTERTERMS)
            return dSigma;

        if (region.hasSoft)
        {
            const double softCT = softCounterterm(point, realChannel, region);
            const double boundaryTerm = softCT * std::log(point.xiMax) * point.xiTilde;

            dSigma += -softCT + boundaryTerm;
        }

        if (region.hasCollinearLeg1)
        {
            const double leg1CollinearCT = leg1CollinearCounterterm(point, realChannel, region, muF2);
            dSigma += -leg1CollinearCT;

            if (region.hasSoft)
            {
                const double leg1SoftCollinearCT = leg1SoftCollinearCounterterm(point, realChannel, region);
                const double boundaryTerm = -leg1SoftCollinearCT * std::log(point.xiMaxLeg1) * point.xiTilde;
                dSigma += leg1SoftCollinearCT + boundaryTerm;
            }
        }

        if (region.hasCollinearLeg2)
        {
            const double leg2CollinearCT = leg2CollinearCounterterm(point, realChannel, region, muF2);
            dSigma += -leg2CollinearCT;

            if (region.hasSoft)
            {
                const double leg2SoftCollinearCT = leg2SoftCollinearCounterterm(point, realChannel, region);
                const double boundaryTerm = -leg2SoftCollinearCT * std::log(point.xiMaxLeg2) * point.xiTilde;
                dSigma += leg2SoftCollinearCT + boundaryTerm;
            }
        }

        return dSigma;
    }

    double RealMinusCTHandler::dSigmaReal(
        const BBarIntegrationPoint& point,
        const RealChannel& realChannel,
        const double muF2
    ) const
    {
        const double amp2 = m_process.realAmp2(point.real, realChannel, point.alphaS);
        const double fluxFactor = 1 / (2.0 * point.real.sHatReal);

        const double f1 = m_config.PDF->xfxQ2(realChannel.id1, point.real.x1, muF2) / point.real.x1;
        const double f2 = m_config.PDF->xfxQ2(realChannel.id2, point.real.x2, muF2) / point.real.x2;
        const double luminosity = f1 * f2;

        return point.jacobian * luminosity * fluxFactor * amp2;
    }

    double RealMinusCTHandler::softCounterterm(
        const BBarIntegrationPoint& point,
        const RealChannel& realChannel,
        const FKSRegion& region
    ) const
    {
        const int emittedID = realChannel.outIDs[region.fksPartonIdx];
        
        if (emittedID != 21)
            throw std::runtime_error("Only Soft gluons have soft singularities");

        const FourVector softGluonVec = {
            1.0,
            std::sqrt(1 - point.y * point.y) * std::sin(point.phi),
            std::sqrt(1 - point.y * point.y) * std::cos(point.phi),
            point.y
        };

        // TODO: Generalize from here
        
        const double B12 = m_config.C_F * point.amp2Born / (2.0 * point.born.sHat);

        const double eikonalFactor = 2.0 * B12 
            * dot(point.born.p1Bar, point.born.p2Bar)
            / dot(point.born.p1Bar, softGluonVec)
            / dot(point.born.p2Bar, softGluonVec);

        const double fluxFactor = 1.0 / (2.0 * point.born.sHat);
        const double luminosity = point.f1Born * point.f2Born;

        const double coupling = 4.0 * PI * point.alphaS;
        const double amp2SoftTimesXiSq = 8.0 * coupling * eikonalFactor;

        // this returns rrrs in sigreal.f
        return point.jacobianOverXiSqSoft * luminosity * fluxFactor * amp2SoftTimesXiSq;
    }

    double RealMinusCTHandler::leg1CollinearCounterterm(
        const BBarIntegrationPoint& point,
        const RealChannel& realChannel,
        const FKSRegion& region,
        const double muF2
    ) const
    {
        // TODO: Generalize this channel selection
        const int emittedID = realChannel.outIDs[region.fksPartonIdx];

        const bool qqChannel = realChannel.id1 == point.bornChannel.id1 && emittedID == 21;
        const bool gqChannel = realChannel.id1 == 21 && realChannel.id2 == point.bornChannel.id2 && emittedID == point.bornChannel.id2;

        const double z  = 1.0 - point.xiLeg1;
        const double x1 = point.born.x1Bar / z;

        const double f1 = m_config.PDF->xfxQ2(realChannel.id1, x1, muF2) / x1;
        const double luminosity = f1 * point.f2Born;

        double kernel = 0.0;
        if (qqChannel)
            kernel = oneMinusZTimesPqq(z);
        else if (gqChannel)
            kernel = oneMinusZTimesPgq(z);
        else
            throw std::runtime_error("This channel has no leg 1 collinear singularity");

        const double fluxFactor = 1.0 / (2.0 * point.born.sHat);

        const double coupling = 4.0 * PI * point.alphaS;
        const double singularFactor = 1.0 / point.xiLeg1 / point.xiLeg1 / (1.0 - point.y);

        const double amp2Collinear =
            4.0 / point.born.sHat
            * coupling
            * kernel
            * singularFactor
            * point.amp2Born;

        // this returns rrrp on line 132 of sigreal.f
        return point.jacobianLeg1 * luminosity * fluxFactor * amp2Collinear;
    }

    double RealMinusCTHandler::leg2CollinearCounterterm(
        const BBarIntegrationPoint& point,
        const RealChannel& realChannel,
        const FKSRegion& region,
        const double muF2
    ) const
    {
        // TODO: Generalize this channel selection
        const int emittedID = realChannel.outIDs[region.fksPartonIdx];

        const bool qqChannel = realChannel.id2 == point.bornChannel.id2 && emittedID == 21;
        const bool gqChannel = realChannel.id2 == 21 && realChannel.id1 == point.bornChannel.id1 && emittedID == point.bornChannel.id1;

        const double z  = 1.0 - point.xiLeg2;
        const double x2 = point.born.x2Bar / z;

        const double f2 = m_config.PDF->xfxQ2(realChannel.id2, x2, muF2) / x2;
        const double luminosity = point.f1Born * f2;

        // const FourVector kPerp = { 0.0, std::sin(point.phi), std::cos(point.phi), 0.0 };

        double kernel = 0.0;
        if (qqChannel)
            kernel = oneMinusZTimesPqq(z);
        else if (gqChannel)
            kernel = oneMinusZTimesPgq(z);
        else
            throw std::runtime_error("This channel has no leg 2 collinear singularity");

        const double fluxFactor = 1.0 / (2.0 * point.born.sHat);

        const double coupling = 4.0 * PI * point.alphaS;
        const double singularFactor = 1.0 / point.xiLeg2 / point.xiLeg2 / (1.0 + point.y);

        const double amp2Collinear =
            4.0 / point.born.sHat
            * coupling
            * kernel
            * singularFactor
            * point.amp2Born;

        // this returns rrrm on line 142 of sigreal.f
        return point.jacobianLeg2 * luminosity * fluxFactor * amp2Collinear;
    }

    double RealMinusCTHandler::leg1SoftCollinearCounterterm(
        const BBarIntegrationPoint& point,
        const RealChannel& realChannel,
        const FKSRegion& region
    ) const
    {
        const int emittedID = realChannel.outIDs[region.fksPartonIdx];

        const bool qqChannel = realChannel.id1 == point.bornChannel.id1 && emittedID == 21;
        const bool gqChannel = realChannel.id1 == 21 && realChannel.id2 == point.bornChannel.id2 && emittedID == point.bornChannel.id2;

        double kernel = 0.0;
        if (qqChannel)
            kernel = oneMinusZTimesPqqAtOne();
        else if (gqChannel)
            kernel = oneMinusZTimesPgqAtOne();
        else
            throw std::runtime_error("This channel has no leg 1 soft collinear singularity");

        const double fluxFactor = 1.0 / (2.0 * point.born.sHat);
        const double luminosity = point.f1Born * point.f2Born;

        const double coupling = 4.0 * PI * point.alphaS;
        const double singularFactorTimesXiSq = 1.0 / (1.0 - point.y);

        const double amp2SoftCollinearTimesXiSq =
            4.0 / point.born.sHat
            * coupling
            * kernel
            * singularFactorTimesXiSq
            * point.amp2Born;
            
        // this returns rrrps on line 134 of sigreal.f
        return point.jacobianOverXiSqSoft * luminosity * fluxFactor * amp2SoftCollinearTimesXiSq;
    }

    double RealMinusCTHandler::leg2SoftCollinearCounterterm(
        const BBarIntegrationPoint& point,
        const RealChannel& realChannel,
        const FKSRegion& region
    ) const
    {
        const int emittedID = realChannel.outIDs[region.fksPartonIdx];

        const bool qqChannel = realChannel.id2 == point.bornChannel.id2 && emittedID == 21;
        const bool gqChannel = realChannel.id2 == 21 && realChannel.id1 == point.bornChannel.id1 && emittedID == point.bornChannel.id1;

        // const FourVector kPerp = { 0.0, std::sin(point.phi), std::cos(point.phi), 0.0 };

        double kernel = 0.0;
        if (qqChannel)
            kernel = oneMinusZTimesPqqAtOne();
        else if (gqChannel)
            kernel = oneMinusZTimesPgqAtOne();
        else
            throw std::runtime_error("This channel has no leg 2 soft collinear singularity");

        const double fluxFactor = 1.0 / (2.0 * point.born.sHat);
        const double luminosity = point.f1Born * point.f2Born;

        const double coupling = 4.0 * PI * point.alphaS;
        const double singularFactorTimesXiSq = 1.0 / (1.0 + point.y);

        const double amp2SoftCollinearTimesXiSq =
            4.0 / point.born.sHat
            * coupling
            * kernel
            * singularFactorTimesXiSq
            * point.amp2Born;
            
        // this returns rrrms on line 144 of sigreal.f
        return point.jacobianOverXiSqSoft * luminosity * fluxFactor * amp2SoftCollinearTimesXiSq;
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
   
} // namespace powheg_dy
