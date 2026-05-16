#include "real_minus_ct_handler.h"

#include "powheg_dy/process.h"

namespace powheg_dy
{
    const double GEV_SQ_TO_PB = 0.389379338e9;

    BBarIntegrationPoint RealMinusCTHandler::generateIntegrationPoint(
        const BornPhSpPt& born, 
        const BornChannel& bornChannel, 
        const double amp2Born, 
        const double muF2,
        const double muR2,
        const std::array<double, 3>& unitCube,
        const ISRRealPhaseSpace& phaseSpace   
    ) const
    {
        BBarIntegrationPoint point;

        point.u1 = unitCube[0];
        point.u2 = unitCube[1];
        point.u3 = unitCube[2];

        point.born = born;
        point.bornChannel = bornChannel;

        point.alphaS = m_config.alphaS(muR2);
        point.amp2Born = amp2Born;

        /////////////////////////////// For the Collinear Remnants ///////////////////////////////

        point.zLeg1CollRemn = born.x1Bar + (1.0 - born.x1Bar) * point.u1;
        point.zLeg2CollRemn = born.x2Bar + (1.0 - born.x2Bar) * point.u1;

        point.f1Born = m_config.PDF->xfxQ2(bornChannel.id1, born.x1Bar, muF2) / born.x1Bar; 
        point.f2Born = m_config.PDF->xfxQ2(bornChannel.id2, born.x2Bar, muF2) / born.x2Bar; 

        point.f1RealQ = m_config.PDF->xfxQ2(bornChannel.id1, born.x1Bar / point.zLeg1CollRemn, muF2) / (born.x1Bar / point.zLeg1CollRemn); 
        point.f2RealQ = m_config.PDF->xfxQ2(bornChannel.id2, born.x2Bar / point.zLeg2CollRemn, muF2) / (born.x2Bar / point.zLeg2CollRemn); 

        point.f1RealG = m_config.PDF->xfxQ2(21, born.x1Bar / point.zLeg1CollRemn, muF2) / (born.x1Bar / point.zLeg1CollRemn); 
        point.f2RealG = m_config.PDF->xfxQ2(21, born.x2Bar / point.zLeg2CollRemn, muF2) / (born.x2Bar / point.zLeg2CollRemn); 
        
        ////////////////////////////////// For the Counterterms //////////////////////////////////
        
        const double TINY_XI = 0.0e-6;
        const double TINY_Y = 0.0e-6;

        double jacobian = 1.0;

        point.xiTilde = (3.0 - 2.0 * point.u1) * point.u1 * point.u1;
        point.xiTilde = point.xiTilde * (1.0 - 2.0 * TINY_XI) + TINY_XI;
        jacobian *= 6.0 * (1.0 - point.u1) * point.u1;

        point.y = 1.0 - 2.0 * point.u2;
        jacobian *= 2.0;
        jacobian *= 1.5 * (1.0 - point.y * point.y);
        point.y = 1.5 * (point.y - std::pow(point.y, 3) / 3) * (1.0 - TINY_Y);

        point.phi = 2.0 * PI * point.u3;
        jacobian *= 2.0 * PI;

        point.xiMax = phaseSpace.xiMax(point.born, point.y);
        point.xiMaxLeg1 = 1.0 - point.born.x1Bar;
        point.xiMaxLeg2 = 1.0 - point.born.x2Bar;

        point.xi = point.xiTilde * point.xiMax;
        point.xiLeg1 = point.xiTilde * point.xiMaxLeg1;
        point.xiLeg2 = point.xiTilde * point.xiMaxLeg2;

        point.real = phaseSpace.reconstruct(born, { point.xi, point.y, point.phi });

        point.jacobianOverXi = jacobian * point.real.radJacobian / point.xi;
        point.jacobianOverXiSoft = jacobian * point.born.sHat / std::pow(4.0 * PI, 3);
        point.jacobianOverXiLeg1 = jacobian * (point.born.sHat / (1.0 - point.xiLeg1)) / std::pow(4.0 * PI, 3) / (1.0 - point.xiLeg1);
        point.jacobianOverXiLeg2 = jacobian * (point.born.sHat / (1.0 - point.xiLeg2)) / std::pow(4.0 * PI, 3) / (1.0 - point.xiLeg2);

        return point;
    }

    double RealMinusCTHandler::dSigmaRealMinusCT(
        const BBarIntegrationPoint& point,
        const RealChannel& realChannel,
        const FKSRegion& region,
        const double muF2,
        const double muR2,
        const ISRRealPhaseSpace& phaseSpace
    ) const
    {
        if (region.hasCollinearFinal)
            throw std::runtime_error("FSR not yet implemented.");

        double dSigma = 0.0;
        
        if (m_config.BTILDE_USE_REAL)
            dSigma += dSigmaReal(point, realChannel, muF2, muR2);
        
        if (!m_config.BTILDE_USE_COUNTERTERMS)
            return dSigma;

        if (region.hasSoft)
        {
            if (false)
            {
                BBarIntegrationPoint pointCp = point;

                for (double n = 1.0e2; n <= 1.0e8; n*=10)
                {
                    auto softPoint = generateIntegrationPoint(point.born, point.bornChannel, point.amp2Born, muF2, muR2, {1.0 / n, point.u2, point.u3 }, phaseSpace);
                    const double a = dSigmaReal(softPoint, realChannel, muF2, muR2);
                    const double b = softCounterterm(softPoint, realChannel, region);
                    if (b == 0) assert(a == 0); else Log::info << a / b << ", a = " << a * GEV_SQ_TO_PB << ", b = " << b * GEV_SQ_TO_PB << Log::endl;
                }
            }
            else
            {
                const auto softPoint = generateIntegrationPoint(point.born, point.bornChannel, point.amp2Born, muF2, muR2, { 10e-6, point.u2, point.u3 }, phaseSpace);
                const double a = dSigmaReal(softPoint, realChannel, muF2, muR2);
                const double b = softCounterterm(softPoint, realChannel, region);
                if (b == 0) assert(a == 0); else if(!(std::abs(a / b - 1.0) < 2e-2)) Log::err << "Soft limit ratio is " << a / b << Log::endl;
            }

            const double softCT = softCounterterm(point, realChannel, region);
            const double boundaryTerm = softCT * std::log(point.xiMax) * point.xiTilde;

            dSigma += -softCT + boundaryTerm;
        }

        if (region.hasCollinearLeg1)
        {
            if (false)
            {
                for (double n = 1.0e1; n <= 1.0e5; n*=10)
                {
                    auto leg1CollinearPoint = generateIntegrationPoint(point.born, point.bornChannel, point.amp2Born, muF2, muR2, {0.25, 1.0 / n, point.u3 }, phaseSpace);
                    const double a = dSigmaReal(leg1CollinearPoint, realChannel, muF2, muR2);
                    const double b = leg1CollinearCounterterm(leg1CollinearPoint, realChannel, region, muF2);
                    if (b == 0) assert(a == 0); else Log::info << a / b << ", a = " << a * GEV_SQ_TO_PB << ", b = " << b * GEV_SQ_TO_PB << Log::endl;
                }
            }
            else
            {
                const auto leg1CollinearPoint = generateIntegrationPoint(point.born, point.bornChannel, point.amp2Born, muF2, muR2, {0.25, 10e-6, point.u3 }, phaseSpace);
                const double a = dSigmaReal(leg1CollinearPoint, realChannel, muF2, muR2);
                const double b = leg1CollinearCounterterm(leg1CollinearPoint, realChannel, region, muF2);
                if (b == 0) assert(a == 0); else if(!(std::abs(a / b - 1.0) < 2e-2)) Log::err << "Leg 1 collinear limit ratio is " << a / b << Log::endl;
            }

            const double leg1CollinearCT = leg1CollinearCounterterm(point, realChannel, region, muF2);
            dSigma += -leg1CollinearCT;

            if (region.hasSoft)
            {
                if (false)
                {
                    for (double n = 1.0e1; n <= 1.0e5; n*=10)
                    {
                        auto leg1SoftCollinearPoint = generateIntegrationPoint(point.born, point.bornChannel, point.amp2Born, muF2, muR2, {1.0 / n, 1.0 / n, point.u3 }, phaseSpace);
                        const double a = dSigmaReal(leg1SoftCollinearPoint, realChannel, muF2, muR2);
                        const double b = leg1SoftCollinearCounterterm(leg1SoftCollinearPoint, realChannel, region, muF2);
                        if (b == 0) assert(a == 0); else Log::info << a / b << ", a = " << a * GEV_SQ_TO_PB << ", b = " << b * GEV_SQ_TO_PB << Log::endl;
                    }
                }
                else
                {
                    const auto leg1SoftCollinearPoint = generateIntegrationPoint(point.born, point.bornChannel, point.amp2Born, muF2, muR2, {10e-6, 10e-6, point.u3 }, phaseSpace);
                    const double a = dSigmaReal(leg1SoftCollinearPoint, realChannel, muF2, muR2);
                    const double b = leg1SoftCollinearCounterterm(leg1SoftCollinearPoint, realChannel, region, muF2);
                    if (b == 0) assert(a == 0); else if(!(std::abs(a / b - 1.0) < 2e-2)) Log::err << "Leg 1 soft collinear limit ratio is " << a / b << Log::endl;
                }

                const double leg1SoftCollinearCT = leg1SoftCollinearCounterterm(point, realChannel, region, muF2);
                const double boundaryTerm = -leg1SoftCollinearCT * std::log(point.xiMaxLeg1) * point.xiTilde;
                dSigma += leg1SoftCollinearCT + boundaryTerm;
            }
        }

        if (region.hasCollinearLeg2)
        {
            if (false)
            {
                BBarIntegrationPoint pointCp = point;

                for (double n = 1.0e1; n <= 1.0e5; n*=10)
                {
                    auto leg2CollinearPoint = generateIntegrationPoint(point.born, point.bornChannel, point.amp2Born, muF2, muR2, {0.25, 1.0 - 1.0 / n, point.u3 }, phaseSpace);
                    const double a = dSigmaReal(leg2CollinearPoint, realChannel, muF2, muR2);
                    const double b = leg2CollinearCounterterm(leg2CollinearPoint, realChannel, region, muF2);
                    if (b == 0) assert(a == 0); else Log::info << a / b << ", a = " << a * GEV_SQ_TO_PB << ", b = " << b * GEV_SQ_TO_PB << Log::endl;
                }
            }
            else
            {
                const auto leg2CollinearPoint = generateIntegrationPoint(point.born, point.bornChannel, point.amp2Born, muF2, muR2, {0.25, 1.0 - 10e-6, point.u3 }, phaseSpace);
                const double a = dSigmaReal(leg2CollinearPoint, realChannel, muF2, muR2);
                const double b = leg2CollinearCounterterm(leg2CollinearPoint, realChannel, region, muF2);
                if (b == 0) assert(a == 0); else if(!(std::abs(a / b - 1.0) < 2e-2)) Log::err << "Leg 2 collinear limit ratio is " << a / b << Log::endl;
            }

            const double leg2CollinearCT = leg2CollinearCounterterm(point, realChannel, region, muF2);
            dSigma += -leg2CollinearCT;

            if (region.hasSoft)
            {
                if (false)
                {
                    for (double n = 1.0e1; n <= 1.0e5; n*=10)
                    {
                        auto leg2SoftCollinearPoint = generateIntegrationPoint(point.born, point.bornChannel, point.amp2Born, muF2, muR2, {1.0 / n, 1.0 - 1.0 / n, point.u3 }, phaseSpace);
                        const double a = dSigmaReal(leg2SoftCollinearPoint, realChannel, muF2, muR2);
                        const double b = leg2SoftCollinearCounterterm(leg2SoftCollinearPoint, realChannel, region, muF2);
                        if (b == 0) assert(a == 0); else Log::info << a / b << ", a = " << a * GEV_SQ_TO_PB << ", b = " << b * GEV_SQ_TO_PB << Log::endl;
                    }
                }
                else
                {
                    const auto leg2SoftCollinearPoint = generateIntegrationPoint(point.born, point.bornChannel, point.amp2Born, muF2, muR2, {10e-6, 1.0 - 10e-6, point.u3 }, phaseSpace);
                    const double a = dSigmaReal(leg2SoftCollinearPoint, realChannel, muF2, muR2);
                    const double b = leg2SoftCollinearCounterterm(leg2SoftCollinearPoint, realChannel, region, muF2);
                    if (b == 0) assert(a == 0); else if(!(std::abs(a / b - 1.0) < 2e-2)) Log::err << "Leg 2 soft collinear limit ratio is " << a / b << Log::endl;
                }

                const double leg2SoftCollinearCT = leg2SoftCollinearCounterterm(point, realChannel, region, muF2);
                const double boundaryTerm = -leg2SoftCollinearCT * std::log(point.xiMaxLeg2) * point.xiTilde;
                dSigma += leg2SoftCollinearCT + boundaryTerm;
            }
        }

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

        // This is rr(j) in sigreal.f
        const double rr = (1.0 - point.y * point.y) * point.xi * point.xi * f1 * f2 * fluxFactor * amp2;

        // This return what is rrr in sigreal.f
        return rr * point.born.jacobian * point.jacobianOverXi / (1.0 - point.y * point.y) / point.xiTilde;
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

        double softalrResult = 4.0 * eikonalFactor * (1.0 - point.y * point.y) / point.born.sHat;
        softalrResult *= 4.0 * PI * point.alphaS;

        // this is rs(j) in sigreal.f
        const double softResult = softalrResult * point.f1Born * point.f2Born;

        // TODO: Apply S_ij factors if more than one soft singularity

        // this returns rrrs in sigreal.f
        return softResult * point.born.jacobian * point.jacobianOverXiSoft / (1.0 - point.y * point.y) / point.xiTilde;
        
        // \frac{2(p1\cdot p2)}{(p1\cdot k)(p2\cdot k)} * \xi^2(1 - y^2)
        // const double eikonalScaled =  16.0 / point.born.sHat;
        // const double softAmp2Scaled = 4.0 * PI * point.alphaS * m_config.C_F * eikonalScaled * point.amp2Born;
        // const double fluxFactor = 1.0 / (2.0 * point.born.sHat);
        // const double luminosity = point.f1Born * point.f2Born;

        // const double dSigmaSoftScaled = luminosity * fluxFactor * softAmp2Scaled;

        // const double jacobianOverXi = point.born.jacobian * point.born.sHat / std::pow(4.0 * PI, 3);

        // // This is the counterterm contribution to ADD to dSigmaReal.
        // // TODO: Move unit cube jacobian here!!
        // return - point.unitCubeJacobian * jacobianOverXi * dSigmaSoftScaled / (point.xi * (1.0 - point.y * point.y)); 
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

        // const FourVector kPerp = { 0.0, std::sin(point.phi), std::cos(point.phi), 0.0 };

        double kernel = 0.0;
        if (qqChannel)
            kernel = oneMinusZTimesPqq(z);
        else if (gqChannel)
            kernel = oneMinusZTimesPgq(z);
        else
            throw std::runtime_error("This channel has no leg 1 collinear singularity");

        const double ap = kernel * point.amp2Born / (2.0 * point.born.sHat);
        const double alr = ap / (point.born.p1Bar.e * point.born.p2Bar.e / z) * 2 * 4 * PI * point.alphaS;
        const double rp = alr * f1 * point.f2Born;

        const double rrrp = rp * point.born.jacobian * point.jacobianOverXiLeg1 / (1.0 - point.y) / point.xiTilde / 2.0;

        return rrrp;
        
        // const double sHatColl = point.born.sHat / z;
        // const double splitting = qqChannel ? oneMinusZTimesPqq(z) : gqChannel ? oneMinusZTimesPgq(z) : 0.0;
        // const double amp2CollScaled = 32.0 * PI * point.alphaS * splitting / sHatColl * point.amp2Born;

        // const double fluxFactor = 1.0 / (2.0 * point.born.sHat);
        // const double radJacobianOverXi = point.born.sHat / (std::pow(4.0 * PI, 3) * (1.0 - xiPlus) * (1.0 - xiPlus));
        // const double unitCubeJacobian = 4.0 * PI * (1.0 - point.born.x1Bar);
        // const double jacobianOverXi = unitCubeJacobian * radJacobianOverXi * point.born.jacobian;

        // const double dSigmaCollScaled = luminosity * jacobianOverXi * fluxFactor * amp2CollScaled;

        // // This is the counterterm contribution to ADD to dSigmaReal.
        // return -dSigmaCollScaled / (2.0 * xiPlus * (1.0 - point.y));
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

        const double ap = kernel * point.amp2Born / (2.0 * point.born.sHat);
        const double alr = ap / (point.born.p1Bar.e * point.born.p2Bar.e / z) * 2 * 4 * PI * point.alphaS;
        const double rm = alr * point.f1Born * f2;

        const double rrrm = rm * point.born.jacobian * point.jacobianOverXiLeg2 / (1.0 + point.y) / point.xiTilde / 2.0;

        return rrrm;
        
        // const double sHatColl = point.born.sHat / z;
        // const double splitting = qqChannel ? oneMinusZTimesPqq(z) : gqChannel ? oneMinusZTimesPgq(z) : 0.0;
        // const double amp2CollScaled = 32.0 * PI * point.alphaS * splitting / sHatColl * point.amp2Born;

        // const double fluxFactor = 1.0 / (2.0 * point.born.sHat);
        // const double radJacobianOverXi = point.born.sHat / (std::pow(4.0 * PI, 3) * (1.0 - xiPlus) * (1.0 - xiPlus));
        // const double unitCubeJacobian = 4.0 * PI * (1.0 - point.born.x1Bar);
        // const double jacobianOverXi = unitCubeJacobian * radJacobianOverXi * point.born.jacobian;

        // const double dSigmaCollScaled = luminosity * jacobianOverXi * fluxFactor * amp2CollScaled;

        // // This is the counterterm contribution to ADD to dSigmaReal.
        // return -dSigmaCollScaled / (2.0 * xiPlus * (1.0 - point.y));
    }

    // double RealMinusCTHandler::leg2CollinearCounterterm(
    //     const BBarIntegrationPoint& point,
    //     const RealChannel& realChannel,
    //     const FKSRegion& region,
    //     const double muF2
    // ) const
    // {
    //     const int emittedID = realChannel.outIDs[region.fksPartonIdx];

    //     // enforce the theta function in the integral to avoid x2Bar/z > 1
    //     if (point.xi > 1.0 - point.born.x2Bar)
    //         return 0.0;

    //     const double z  = 1.0 - point.xi;
    //     const double x2 = point.born.x2Bar / z;

    //     // TODO: Generalize from here
        
    //     const bool qqChannel = realChannel.id2 == point.bornChannel.id2 && emittedID == 21;
    //     const bool gqChannel = realChannel.id2 == 21 && realChannel.id1 == point.bornChannel.id1 && emittedID == point.bornChannel.id1;

    //     if (!qqChannel && !gqChannel)
    //         throw std::runtime_error("This channel has no leg 2 collinear singularity");

    //     const double f2 = m_config.PDF->xfxQ2(realChannel.id2, x2, muF2) / x2;
    //     const double luminosity = point.f1Born * f2;

    //     const double sHatColl = point.born.sHat / z;
    //     const double splitting = qqChannel ? oneMinusZTimesPqq(z) : gqChannel ? oneMinusZTimesPgq(z) : 0.0;
    //     const double amp2CollScaled = 32.0 * PI * point.alphaS * splitting / sHatColl * point.amp2Born;

    //     const double fluxFactor = 1.0 / (2.0 * point.born.sHat);
    //     const double jacobian = point.real.radJacobian * point.born.jacobian;

    //     const double dSigmaCollScaled = luminosity * jacobian * fluxFactor * amp2CollScaled;

    //     // This is the counterterm contribution to ADD to dSigmaReal.
    //     return -dSigmaCollScaled;
    // }

    double RealMinusCTHandler::leg1SoftCollinearCounterterm(
        const BBarIntegrationPoint& point,
        const RealChannel& realChannel,
        const FKSRegion& region,
        const double muF2
    ) const
    {
        const int emittedID = realChannel.outIDs[region.fksPartonIdx];

        const bool qqChannel = realChannel.id1 == point.bornChannel.id1 && emittedID == 21;
        const bool gqChannel = realChannel.id1 == 21 && realChannel.id2 == point.bornChannel.id2 && emittedID == point.bornChannel.id2;

        // const FourVector kPerp = { 0.0, std::sin(point.phi), std::cos(point.phi), 0.0 };

        double kernel = 0.0;
        if (qqChannel)
            kernel = oneMinusZTimesPqqAtOne();
        else if (gqChannel)
            kernel = oneMinusZTimesPgqAtOne();
        else
            throw std::runtime_error("This channel has no leg 1 soft collinear singularity");

        const double ap = kernel * point.amp2Born / (2.0 * point.born.sHat);
        const double alr = ap / (point.born.p1Bar.e * point.born.p2Bar.e) * 2 * 4 * PI * point.alphaS;
        const double rps = alr * point.f1Born * point.f2Born;

        const double rrrps = rps * point.born.jacobian * point.jacobianOverXiSoft / (1.0 - point.y) / point.xiTilde / 2.0;

        return rrrps;

        // const int emittedID = realChannel.outIDs[region.fksPartonIdx];

        // // TODO: Generalize from here

        // const bool qqChannel = realChannel.id1 == point.bornChannel.id1 && emittedID == 21;

        // if (!qqChannel)
        //     throw std::runtime_error("This channel has no leg 1 collinear singularity");

        // const double xiPlus = point.u1 * (1.0 - point.born.x1Bar);

        // const double luminosity = point.f1Born * point.f2Born;

        // const double splitting = qqChannel ? oneMinusZTimesPqqAtOne() : 0.0;
        // const double amp2SoftCollScaled = 32.0 * PI * point.alphaS * splitting / point.born.sHat * point.amp2Born;

        // const double fluxFactor = 1.0 / (2.0 * point.born.sHat);

        // const double unitCubeJacobian = 4.0 * PI * (1.0 - point.born.x1Bar);
        // const double radJacobianOverXi = point.born.sHat / std::pow(4.0 * PI, 3);
        // const double jacobianOverXi = unitCubeJacobian * radJacobianOverXi * point.born.jacobian;

        // const double dSigmaSoftCollScaled = luminosity * jacobianOverXi * fluxFactor * amp2SoftCollScaled;

        // // Soft-collinear overlap is added back.
        // return dSigmaSoftCollScaled / (2.0 * xiPlus * (1.0 - point.y));
    }

    double RealMinusCTHandler::leg2SoftCollinearCounterterm(
        const BBarIntegrationPoint& point,
        const RealChannel& realChannel,
        const FKSRegion& region,
        const double muF2
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

        const double ap = kernel * point.amp2Born / (2.0 * point.born.sHat);
        const double alr = ap / (point.born.p1Bar.e * point.born.p2Bar.e) * 2 * 4 * PI * point.alphaS;
        const double rps = alr * point.f1Born * point.f2Born;

        const double rrrps = rps * point.born.jacobian * point.jacobianOverXiSoft / (1.0 + point.y) / point.xiTilde / 2.0;

        return rrrps;
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
