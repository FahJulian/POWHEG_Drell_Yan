#include "real_minus_ct_handler.h"

#include "powheg_dy/process.h"

namespace powheg_dy
{
    double RealMinusCTHandler::dSigmaRealMinusCT(
        const BBarIntegrationPoint& point,
        const RealChannel& realChannel,
        const FKSRegion& region,
        const double muF2,
        const double muR2,
        const ISRRealPhaseSpace& phaseSpace
    ) const
    {
        double dSigma = dSigmaReal(point, realChannel, muF2, muR2);

        if (region.hasCollinearFinal)
            throw std::runtime_error("FSR not yet implemented.");

        if (region.hasSoft)
        {
            if (true)
            {
                BBarIntegrationPoint pointCp = point;

                for (double n = 1.0e4; n <= 1.0e12; n*=10)
                {
                    pointCp.u1 = 1.0 / n;

                    pointCp.y = -1.0 + 2.0 * pointCp.u2;
                    pointCp.xi = phaseSpace.xiMax(pointCp.born, pointCp.y) * pointCp.u1;
                    pointCp.phi = 2.0 * PI * pointCp.u3;

                    pointCp.real = phaseSpace.reconstruct(pointCp.born, { pointCp.xi, pointCp.y, pointCp.phi });

                    Log::info << dSigmaReal(pointCp, realChannel, muF2, muR2) * (pointCp.xi * pointCp.xi * (1.0 - pointCp.y * pointCp.y))
                            / -softCTScaled(pointCp, realChannel, region) 
                        << Log::endl;
                }
            }

            const double softScaled = softCTScaled(point, realChannel, region);
            dSigma += softScaled / (point.xi * point.xi * (1.0 - point.y * point.y));
        }

        if (region.hasCollinearLeg1)
        {
            if (false)
            {
                BBarIntegrationPoint pointCp = point;

                for (double n = 1.0e4; n <= 1.0e10; n*=10)
                {
                    pointCp.y = 1.0 - 1.0 / n;
                    pointCp.xi = phaseSpace.xiMax(pointCp.born, pointCp.y) * pointCp.u1;
                    pointCp.phi = 2.0 * PI * pointCp.u3;

                    pointCp.real = phaseSpace.reconstruct(pointCp.born, { pointCp.xi, pointCp.y, pointCp.phi });

                    Log::info << dSigmaReal(pointCp, realChannel, muF2, muR2) * (pointCp.xi * pointCp.xi * (1.0 - pointCp.y * pointCp.y))
                            / -leg1CollinearCTScaled(pointCp, realChannel, region, muF2) 
                        << Log::endl;
                }
            }

            const double leg1CollinearScaled = leg1CollinearCTScaled(point, realChannel, region, muF2);
            dSigma += leg1CollinearScaled / (2.0 * point.xi * point.xi * (1.0 - point.y));

            if (region.hasSoft)
            {
                if (false)
                {
                    BBarIntegrationPoint pointCp = point;

                    for (double n = 1.0e4; n <= 1.0e10; n*=10)
                    {
                        pointCp.y = 1.0 - 1.0 / n;
                        pointCp.xi = 1.0 / n;
                        pointCp.phi = 2.0 * PI * pointCp.u3;

                        pointCp.real = phaseSpace.reconstruct(pointCp.born, { pointCp.xi, pointCp.y, pointCp.phi });

                        Log::info << dSigmaReal(pointCp, realChannel, muF2, muR2) * (pointCp.xi * pointCp.xi * (1.0 - pointCp.y * pointCp.y))
                                / leg1SoftCollinearCTScaled(pointCp, realChannel, region, muF2) 
                            << Log::endl;
                    }
                }

                const double leg1SoftCollinearScaled = leg1SoftCollinearCTScaled(point, realChannel, region, muF2);
                dSigma += leg1SoftCollinearScaled / (2.0 * point.xi * point.xi * (1.0 - point.y));
            }
        }

        if (region.hasCollinearLeg2)
        {
            if (false)
            {
                BBarIntegrationPoint pointCp = point;

                for (double n = 1.0e4; n <= 1.0e10; n*=10)
                {
                    pointCp.y = -1.0 + 1.0 / n;
                    pointCp.xi = phaseSpace.xiMax(pointCp.born, pointCp.y) * pointCp.u1;
                    pointCp.phi = 2.0 * PI * pointCp.u3;

                    pointCp.real = phaseSpace.reconstruct(pointCp.born, { pointCp.xi, pointCp.y, pointCp.phi });

                    Log::info << dSigmaReal(pointCp, realChannel, muF2, muR2) * (pointCp.xi * pointCp.xi * (1.0 - pointCp.y * pointCp.y))
                            / -leg2CollinearCTScaled(pointCp, realChannel, region, muF2) 
                        << Log::endl;
                }
            }

            const double leg2CollinearScaled = leg2CollinearCTScaled(point, realChannel, region, muF2);
            dSigma += leg2CollinearScaled / (2.0 * point.xi * point.xi * (1.0 + point.y));

            if (region.hasSoft)
            {
                if (false)
                {
                    BBarIntegrationPoint pointCp = point;

                    for (double n = 1.0e4; n <= 1.0e10; n*=10)
                    {
                        pointCp.y = -1.0 + 1.0 / n;
                        pointCp.xi = 1.0 / n;
                        pointCp.phi = 2.0 * PI * pointCp.u3;

                        pointCp.real = phaseSpace.reconstruct(pointCp.born, { pointCp.xi, pointCp.y, pointCp.phi });

                        Log::info << dSigmaReal(pointCp, realChannel, muF2, muR2) * (pointCp.xi * pointCp.xi * (1.0 - pointCp.y * pointCp.y))
                                / leg2SoftCollinearCTScaled(pointCp, realChannel, region, muF2) 
                            << Log::endl;
                    }
                }

                const double leg2SoftCollinearScaled = leg2SoftCollinearCTScaled(point, realChannel, region, muF2);
                dSigma += leg2SoftCollinearScaled / (2.0 * point.xi * point.xi * (1.0 + point.y));
            }
        }

        // TODO: FSR Collinear 

        dSigma += integralBoundaryRemnant(point, realChannel, region, muF2, phaseSpace);

        return dSigma;
    }

    double RealMinusCTHandler::dSigmaReal(
        const BBarIntegrationPoint& point,
        const RealChannel& realChannel,
        const double muF2,
        const double muR2
    ) const
    {
        const double amp2 = m_process.realAmp2(point.real, point.bornChannel, realChannel, point.alphaS);
        const double fluxFactor = 1 / (2.0 * point.real.sHatReal);

        const double f1 = m_config.PDF->xfxQ2(realChannel.id1, point.real.x1, muF2) / point.real.x1;
        const double f2 = m_config.PDF->xfxQ2(realChannel.id2, point.real.x2, muF2) / point.real.x2;

        const double dSigmaPartonic = point.real.radJacobian * point.born.jacobian * fluxFactor * amp2; 
        const double dSigma = f1 * f2 * dSigmaPartonic;

        return dSigma;
    }

    double RealMinusCTHandler::softCTScaled(
        const BBarIntegrationPoint& point,
        const RealChannel& realChannel,
        const FKSRegion& region
    ) const
    {
        const int emittedID = realChannel.outIDs[region.fksPartonIdx];
        
        if (emittedID != 21)
            throw std::runtime_error("Only Soft gluons have soft singularities");

        // TODO: Generalize from here
        
        // \frac{2(p1\cdot p2)}{(p1\cdot k)(p2\cdot k)} * \xi^2(1 - y^2)
        const double eikonalScaled =  16.0 / point.real.sHatReal;
        const double softAmp2Scaled = 4.0 * PI * point.alphaS * m_config.C_F * eikonalScaled * point.amp2Born;
        const double fluxFactor = 1.0 / (2.0 * point.real.sHatReal);
        const double luminosity = point.f1Born * point.f2Born;
        const double jacobian = point.born.jacobian * point.real.radJacobian;

        const double dSigmaSoftScaled = luminosity * fluxFactor * jacobian * softAmp2Scaled;

        // This is the counterterm contribution to ADD to dSigmaReal.
        return -dSigmaSoftScaled; 
    }

    double RealMinusCTHandler::leg1CollinearCTScaled(
        const BBarIntegrationPoint& point,
        const RealChannel& realChannel,
        const FKSRegion& region,
        const double muF2
    ) const
    {
        const int emittedID = realChannel.outIDs[region.fksPartonIdx];

        // enforce the theta function in the integral to avoid x1Bar/z > 1
        if (point.xi > 1.0 - point.born.x1Bar)
            return 0.0;

        const double z  = 1.0 - point.xi;
        const double x1 = point.born.x1Bar / z;

        // TODO: Generalize from here

        const bool qqChannel = realChannel.id1 == point.bornChannel.id1 && emittedID == 21;
        const bool gqChannel = realChannel.id1 == 21 && realChannel.id2 == point.bornChannel.id2 && emittedID == point.bornChannel.id2;

        if (!qqChannel && !gqChannel)
            throw std::runtime_error("This channel has no leg 1 collinear singularity");

        const double f1 = m_config.PDF->xfxQ2(realChannel.id1, x1, muF2) / x1;
        const double luminosity = f1 * point.f2Born;

        const double sHatColl = point.born.sHat / z;
        const double splitting = qqChannel ? oneMinusZTimesPqq(z) : gqChannel ? oneMinusZTimesPgq(z) : 0.0;
        const double amp2CollScaled = 32.0 * PI * point.alphaS * splitting / sHatColl * point.amp2Born;

        const double fluxFactor = 1.0 / (2.0 * point.born.sHat);
        const double jacobian = point.real.radJacobian * point.born.jacobian;

        const double dSigmaCollScaled = luminosity * jacobian * fluxFactor * amp2CollScaled;

        // This is the counterterm contribution to ADD to dSigmaReal.
        return -dSigmaCollScaled;
    }

    double RealMinusCTHandler::leg2CollinearCTScaled(
        const BBarIntegrationPoint& point,
        const RealChannel& realChannel,
        const FKSRegion& region,
        const double muF2
    ) const
    {
        const int emittedID = realChannel.outIDs[region.fksPartonIdx];

        // enforce the theta function in the integral to avoid x2Bar/z > 1
        if (point.xi > 1.0 - point.born.x2Bar)
            return 0.0;

        const double z  = 1.0 - point.xi;
        const double x2 = point.born.x2Bar / z;

        // TODO: Generalize from here
        
        const bool qqChannel = realChannel.id2 == point.bornChannel.id2 && emittedID == 21;
        const bool gqChannel = realChannel.id2 == 21 && realChannel.id1 == point.bornChannel.id1 && emittedID == point.bornChannel.id1;

        if (!qqChannel && !gqChannel)
            throw std::runtime_error("This channel has no leg 2 collinear singularity");

        const double f2 = m_config.PDF->xfxQ2(realChannel.id2, x2, muF2) / x2;
        const double luminosity = point.f1Born * f2;

        const double sHatColl = point.born.sHat / z;
        const double splitting = qqChannel ? oneMinusZTimesPqq(z) : gqChannel ? oneMinusZTimesPgq(z) : 0.0;
        const double amp2CollScaled = 32.0 * PI * point.alphaS * splitting / sHatColl * point.amp2Born;

        const double fluxFactor = 1.0 / (2.0 * point.born.sHat);
        const double jacobian = point.real.radJacobian * point.born.jacobian;

        const double dSigmaCollScaled = luminosity * jacobian * fluxFactor * amp2CollScaled;

        // This is the counterterm contribution to ADD to dSigmaReal.
        return -dSigmaCollScaled;
    }

    double RealMinusCTHandler::leg1SoftCollinearCTScaled(
        const BBarIntegrationPoint& point,
        const RealChannel& realChannel,
        const FKSRegion& region,
        const double muF2
    ) const
    {
        const int emittedID = realChannel.outIDs[region.fksPartonIdx];

        // TODO: Generalize from here

        const bool qqChannel = realChannel.id1 == point.bornChannel.id1 && emittedID == 21;

        if (!qqChannel)
            throw std::runtime_error("This channel has no leg 1 collinear singularity");

        const double luminosity = point.f1Born * point.f2Born;

        const double splitting = qqChannel ? oneMinusZTimesPqqAtOne() : 0.0;
        const double amp2SoftCollScaled = 32.0 * PI * point.alphaS * splitting / point.born.sHat * point.amp2Born;

        const double fluxFactor = 1.0 / (2.0 * point.born.sHat);
        const double jacobian = point.real.radJacobian * point.born.jacobian;

        const double dSigmaSoftCollScaled = luminosity * jacobian * fluxFactor * amp2SoftCollScaled;

        // This is the counterterm contribution to ADD to dSigmaReal.
        return dSigmaSoftCollScaled;
    }

    double RealMinusCTHandler::leg2SoftCollinearCTScaled(
        const BBarIntegrationPoint& point,
        const RealChannel& realChannel,
        const FKSRegion& region,
        const double muF2
    ) const
    {
        const int emittedID = realChannel.outIDs[region.fksPartonIdx];

        // TODO: Generalize from here

        const bool qqChannel = realChannel.id2 == point.bornChannel.id2 && emittedID == 21;

        if (!qqChannel)
            throw std::runtime_error("This channel has no leg 2 soft collinear singularity");

        const double luminosity = point.f1Born * point.f2Born;

        const double splitting = qqChannel ? oneMinusZTimesPqqAtOne() : 0.0;
        const double amp2SoftCollScaled = 32.0 * PI * point.alphaS * splitting / point.born.sHat * point.amp2Born;

        const double fluxFactor = 1.0 / (2.0 * point.born.sHat);
        const double jacobian = point.real.radJacobian * point.born.jacobian;

        const double dSigmaSoftCollScaled = luminosity * jacobian * fluxFactor * amp2SoftCollScaled;

        // This is the counterterm contribution to ADD to dSigmaReal.
        return dSigmaSoftCollScaled;
    }

    double RealMinusCTHandler::integralBoundaryRemnant(
        const BBarIntegrationPoint& point,
        const RealChannel& realChannel,
        const FKSRegion& region,
        const double muF2,
        const ISRRealPhaseSpace& phaseSpace
    ) const
    {
        if (!region.hasSoft)
            return 0.0;

        const double xi = point.xi;
        const double y  = point.y;

        const double xiTilde = point.u1; // since xi = xiMax(y) * u1

        // If you use the maximal plus prescription, xiCut = 1.
        // If you later introduce a separate xi_cut parameter, put it here.
        const double xiCut = 1.0;

        const double xiMaxFull =
            phaseSpace.xiMax(point.born, y);

        const double xiMaxLeg1 =
            1.0 - point.born.x1Bar;

        const double xiMaxLeg2 =
            1.0 - point.born.x2Bar;

        double remnant = 0.0;

        // Positive soft approximation density.
        //
        // softCTScaled is negative, so use the minus sign.
        const double softDensity =
            -softCTScaled(point, realChannel, region)
            / (xi * xi * (1.0 - y * y));

        remnant += softDensity
            * std::log(xiMaxFull / xiCut)
            * xiTilde;

        if (region.hasCollinearLeg1)
        {
            const double softCollLeg1Density =
                leg1SoftCollinearCTScaled(point, realChannel, region, muF2)
                / (2.0 * xi * xi * (1.0 - y));

            remnant -= softCollLeg1Density
                * std::log(xiMaxLeg1 / xiCut)
                * xiTilde;
        }

        if (region.hasCollinearLeg2)
        {
            const double softCollLeg2Density =
                leg2SoftCollinearCTScaled(point, realChannel, region, muF2)
                / (2.0 * xi * xi * (1.0 + y));

            remnant -= softCollLeg2Density
                * std::log(xiMaxLeg2 / xiCut)
                * xiTilde;
        }

        return remnant;
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
