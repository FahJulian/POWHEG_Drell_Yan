#include "bbar_integrator.h"

#include "powheg_dy/process.h"
#include "powheg_dy/math/rand.h"
#include "powheg_dy/bbar/singular_region.h"
#include "powheg_dy/bbar/collinear_remnant_channel.h"

namespace powheg_dy
{
namespace 
{
    static double p1 = 0.0;
    static double p2 = 0.0;
    static double p3 = 0.0;
    static int n = 0;
    static constexpr double GEV_SQ_TO_PB = 0.389379338e9;
    static constexpr double SECURITY_FACTOR = 2.5;
    static constexpr int MAX_TRIALS = 10000;

    double sPlus(double y)
    {
        return 0.5 * (1.0 + y);
    }

    double sMinus(double y) 
    {
        return 0.5 * (1.0 - y);
    }

    enum RealChannelOld : int
    {
        QQBAR = 0,
        GLUON_LEG1 = 1,
        GLUON_LEG2 = 2,
    };

} // anonymous namespace 
    
    void BBarIntegrator::computeWeightAndSampleChannel(BornPhSpPt& born) const
    {
        const std::array<double, 3> unitCube = { rand(), rand(), rand() };        // For the radiation variables
        // const std::array<double, 4> unitX = { rand(), rand(), rand(), rand() };  // For the collinear remnants

        std::vector<BornChannel> channels = m_process.bornChannels();
        
        std::vector<std::pair<BornPhSpPt, double>> weights = {};
        weights.reserve(channels.size());
            
        double totalAbsWeight = 0.0;
        for (const auto& channel : channels)
        {
            BornPhSpPt bornCopy = born;
            bornCopy.channel = channel;
            m_bornPhaseSpace->reconstructMomenta(bornCopy);

            // const double weight = bTildeOld(bornCopy, unitCube, unitX);
            const double weight = bTilde(bornCopy, unitCube);

            totalAbsWeight += std::abs(weight);

            weights.push_back({ bornCopy, weight });
        }

        // Sample the parton channel by their relative contribution to dSigma
        double u = rand(0.0, totalAbsWeight);
        for (auto& [bornCopy, weight] : weights)
        {
            if (u < std::abs(weight))
            {   
                bornCopy.weight = totalAbsWeight;
                bornCopy.sign = weight >= 0 ? 1 : -1;
                born = bornCopy;

                break;
            }

            u -= std::abs(weight);
        }
    }

    BornPhSpPt BBarIntegrator::sampleAccordingtoBTilde()
    {
        for (int trials = 1; trials <= MAX_TRIALS; trials++)
        {
            double rands[3] = { rand(), rand(), rand() };
            BornPhSpPt born = m_bornPhaseSpace->samplePoint(rands);
            
            computeWeightAndSampleChannel(born);

            // powheg_assert(born.weight <= m_maxWeight, "Born weight " << born.weight << " exceeds max weight " << m_maxWeight);
            if (born.weight > m_maxWeight)
            {
                Log::warn << "Born weight " << born.weight << " exceeds max weight " << m_maxWeight << std::endl;
                m_nEventTrials += trials;
                m_sumSigns += born.sign;
                return born;
            }

            double u = rand();
            if (u < born.weight / m_maxWeight)
            {
                m_nEventTrials += trials;
                m_sumSigns += born.sign;
                return born;
            }
        }

        powheg_assert(false, "Bbar sampling stuck in loop.");
        return { };
    }

    RadiationVariables BBarIntegrator::sampleUniformRad(const BornPhSpPt& born, const std::array<double, 3>& unitCube) const
    {
        const double y = -1.0 + 2.0 * unitCube[1];

        const double xiMax = m_realPhaseSpace->xiMax(born, y);
        const double xi = xiMax * unitCube[0];

        const double phi = 2.0 * PI * unitCube[2];
        
        return { xi, y, phi, xiMax };
    }

    double BBarIntegrator::unitCubeJacobian(const RadiationVariables& rad) const
    {
        return 4.0 * PI * rad.xiMax;
    }

    double BBarIntegrator::luminosity(const BornPhSpPt& born, const RadiationVariables& rad, int id1, int id2, double muF2) const
    {
        const double xi = rad.xi;
        const double y  = rad.y;

        const double a = 2.0 - xi * (1.0 - y);
        const double b = 2.0 - xi * (1.0 + y);

        const double sqrtOneMinusXi = std::sqrt(1.0 - xi);

        const double x1 = born.x1Bar / sqrtOneMinusXi * std::sqrt(a / b);
        const double x2 = born.x2Bar / sqrtOneMinusXi * std::sqrt(b / a);

        const double fLeg1 = m_config.PDF->xfxQ2(id1, x1, muF2) / x1; 
        const double fLeg2 = m_config.PDF->xfxQ2(id2, x2, muF2) / x2; 

        return fLeg1 * fLeg2;
    }

    double BBarIntegrator::AqqbarPlus(const RealPhSpPt& real, double muR2) const
    {
        const BornChannel& bornChannel = real.underlyingBorn.channel;
        RealChannel channel = { 
            .id1    = bornChannel.id1, 
            .id2    = bornChannel.id2, 
            .outIDs = { bornChannel.outIDs[0], bornChannel.outIDs[1], 21 } 
        };

        return sPlus(real.rad.y) * (1.0 - real.rad.y) * real.rad.xi * real.rad.xi 
            * m_process.realAmp2(real, channel, m_config.alphaS(muR2)) / real.sHatReal;
    }

    double BBarIntegrator::AqqbarMinus(const RealPhSpPt& real, double muR2) const
    {
        const BornChannel& bornChannel = real.underlyingBorn.channel;
        RealChannel channel = { 
            .id1    = bornChannel.id1, 
            .id2    = bornChannel.id2, 
            .outIDs = { bornChannel.outIDs[0], bornChannel.outIDs[1], 21 } 
        };

        return sMinus(real.rad.y) * (1.0 + real.rad.y) * real.rad.xi * real.rad.xi 
            * m_process.realAmp2(real, channel, m_config.alphaS(muR2)) / real.sHatReal;
    }

    double BBarIntegrator::AqqbarPlusLimitAware(const BornPhSpPt& born, const RadiationVariables& rad, double muR2) const
    {
        const bool soft = (rad.xi == 0.0);
        const bool coll = (rad.y  == 1.0);

        if (soft && coll)
            return AqqbarPlusSoftCollinearLimit(born, muR2);
        else if (soft)
            return AqqbarPlusSoftLimit(born, rad, muR2);
        else if (coll)
            return AqqbarPlusCollinearLimit(born, rad, muR2);
        else
        {
            RealPhSpPt real = m_realPhaseSpace->reconstruct(born, rad);
            return AqqbarPlus(real, muR2);
        }
    }

    double BBarIntegrator::AqqbarMinusLimitAware(const BornPhSpPt& born, const RadiationVariables& rad, double muR2) const
    {
        const bool soft = (rad.xi ==  0.0);
        const bool coll = (rad.y  == -1.0);

        if (soft && coll)
            return AqqbarMinusSoftCollinearLimit(born, muR2);
        else if (soft)
            return AqqbarMinusSoftLimit(born, rad, muR2);
        else if (coll)
            return AqqbarMinusCollinearLimit(born, rad, muR2);
        else
        {
            RealPhSpPt real = m_realPhaseSpace->reconstruct(born, rad);
            return AqqbarMinus(real, muR2);
        }
    }

    // TODO: Merge the two functions below into one
    // TODO: Create a UnitCube struct with meaningful names
    double BBarIntegrator::RHatqqbarPlus(const BornPhSpPt& born, const std::array<double, 3>& unitCube, double muF2, double muR2) const
    {
        const RadiationVariables rad         = sampleUniformRad(born, unitCube);
        const RadiationVariables radSoft     = sampleUniformRad(born, {0.0,         unitCube[1], unitCube[2]});
        const RadiationVariables radColl     = sampleUniformRad(born, {unitCube[0], 1.0,         unitCube[2]});
        const RadiationVariables radSoftColl = sampleUniformRad(born, {0.0,         1.0,         unitCube[2]});

        const double lum         = luminosity(born, rad,         born.channel.id1, born.channel.id2, muF2);
        const double lumSoft     = luminosity(born, radSoft,     born.channel.id1, born.channel.id2, muF2);
        const double lumColl     = luminosity(born, radColl,     born.channel.id1, born.channel.id2, muF2);
        const double lumSoftColl = luminosity(born, radSoftColl, born.channel.id1, born.channel.id2, muF2);

        const double A_xy = lum         * AqqbarPlusLimitAware(born, rad,         muR2);
        const double A_0y = lumSoft     * AqqbarPlusLimitAware(born, radSoft,     muR2);
        const double A_x1 = lumColl     * AqqbarPlusLimitAware(born, radColl,     muR2);
        const double A_01 = lumSoftColl * AqqbarPlusLimitAware(born, radSoftColl, muR2);

        return (A_xy - A_0y - A_x1 + A_01) / (rad.xi * rad.xi * (1.0 - rad.y));
    }

    double BBarIntegrator::RHatqqbarMinus(const BornPhSpPt& born, const std::array<double, 3>& unitCube, double muF2, double muR2) const
    {
        const RadiationVariables rad         = sampleUniformRad(born, unitCube);
        const RadiationVariables radSoft     = sampleUniformRad(born, {0.0,         unitCube[1], unitCube[2]});
        const RadiationVariables radColl     = sampleUniformRad(born, {unitCube[0], 0.0,         unitCube[2]});
        const RadiationVariables radSoftColl = sampleUniformRad(born, {0.0,         0.0,         unitCube[2]});

        const double lum         = luminosity(born, rad,         born.channel.id1, born.channel.id2, muF2);
        const double lumSoft     = luminosity(born, radSoft,     born.channel.id1, born.channel.id2, muF2);
        const double lumColl     = luminosity(born, radColl,     born.channel.id1, born.channel.id2, muF2);
        const double lumSoftColl = luminosity(born, radSoftColl, born.channel.id1, born.channel.id2, muF2);

        const double A_xy = lum         * AqqbarMinusLimitAware(born, rad,         muR2);
        const double A_0y = lumSoft     * AqqbarMinusLimitAware(born, radSoft,     muR2);
        const double A_xm = lumColl     * AqqbarMinusLimitAware(born, radColl,     muR2);
        const double A_0m = lumSoftColl * AqqbarMinusLimitAware(born, radSoftColl, muR2);

        return (A_xy - A_0y - A_xm + A_0m) / (rad.xi * rad.xi * (1.0 + rad.y));
    }

    double BBarIntegrator::AqqbarPlusCollinearLimit(const BornPhSpPt& born, const RadiationVariables& rad, double muR2) const
    {
        const double z = 1.0 - rad.xi;
        const double prefactor = 16.0 * PI * m_config.C_F / born.sHat / born.sHat;

        return prefactor
            * m_config.alphaS(muR2)
            * m_process.bornAmp2(born)
            * z * (1 + z * z);
    }

    double BBarIntegrator::AqqbarMinusCollinearLimit(const BornPhSpPt& born, const RadiationVariables& rad, double muR2) const
    {
        const double z = 1.0 - rad.xi;
        const double prefactor = 16.0 * PI * m_config.C_F / born.sHat / born.sHat;

        return prefactor
            * m_config.alphaS(muR2)
            * m_process.bornAmp2(born)
            * z * (1 + z * z);
    }

    double BBarIntegrator::AqqbarPlusSoftLimit(const BornPhSpPt& born, const RadiationVariables& rad, double muR2) const
    {
        const double angular = (1.0 - rad.y) * sPlus(rad.y) / (1.0 - rad.y * rad.y);
        const double prefactor = 64.0 * PI * m_config.C_F / born.sHat / born.sHat;

        return prefactor
            * m_config.alphaS(muR2)
            * m_process.bornAmp2(born)
            * angular;
    }

    double BBarIntegrator::AqqbarMinusSoftLimit(const BornPhSpPt& born, const RadiationVariables& rad, double muR2) const
    {
        const double angular = (1.0 + rad.y) * sMinus(rad.y) / (1.0 - rad.y * rad.y);
        const double prefactor = 64.0 * PI * m_config.C_F / born.sHat / born.sHat;

        return prefactor
            * m_config.alphaS(muR2)
            * m_process.bornAmp2(born)
            * angular;
    }

    double BBarIntegrator::AqqbarPlusSoftCollinearLimit(const BornPhSpPt& born, double muR2) const
    {
        const double prefactor = 32.0 * PI * m_config.C_F / born.sHat / born.sHat;
        return prefactor * m_config.alphaS(muR2) * m_process.bornAmp2(born);
    }

    double BBarIntegrator::AqqbarMinusSoftCollinearLimit(const BornPhSpPt& born, double muR2) const
    {
        const double prefactor = 32.0 * PI * m_config.C_F / born.sHat / born.sHat;
        return prefactor * m_config.alphaS(muR2) * m_process.bornAmp2(born);
    }

    double BBarIntegrator::AgluonLeg1Plus(const RealPhSpPt& real, double muR2) const
    {
        const BornChannel& bornChannel = real.underlyingBorn.channel;
        RealChannel channel = { 
            .id1    = 21, 
            .id2    = bornChannel.id2, 
            .outIDs = { bornChannel.outIDs[0], bornChannel.outIDs[1], bornChannel.id2 } 
        };

        return (1.0 - real.rad.y) * m_process.realAmp2(real, channel, m_config.alphaS(muR2)) / real.sHatReal;
    }

    double BBarIntegrator::AgluonLeg2Minus(const RealPhSpPt& real, double muR2) const
    {
        const BornChannel& bornChannel = real.underlyingBorn.channel;
        RealChannel channel = { 
            .id1    = bornChannel.id1, 
            .id2    = 21, 
            .outIDs = { bornChannel.outIDs[0], bornChannel.outIDs[1], bornChannel.id1 } 
        };

        return (1.0 + real.rad.y) * m_process.realAmp2(real, channel, m_config.alphaS(muR2)) / real.sHatReal;
    }

    double BBarIntegrator::AgluonLeg1PlusCollinearLimit(const BornPhSpPt& born, const RadiationVariables& rad, double muR2) const
    {
        const double z = 1.0 - rad.xi;
        const double prefactor = 16.0 * PI * m_config.T_F / born.sHat / born.sHat;

        return m_process.bornAmp2(born)
            * prefactor
            * m_config.alphaS(muR2)
            * z
            * (z * z + (1.0 - z) * (1.0 - z))
            / (1.0 - z);
    }

    double BBarIntegrator::AgluonLeg2MinusCollinearLimit(const BornPhSpPt& born, const RadiationVariables& rad, double muR2) const
    {
        const double z = 1.0 - rad.xi;
        const double prefactor = 16.0 * PI * m_config.T_F / born.sHat / born.sHat;

        return m_process.bornAmp2(born)
            * prefactor
            * m_config.alphaS(muR2)
            * z
            * (z * z + (1.0 - z) * (1.0 - z))
            / (1.0 - z);
    }

    double BBarIntegrator::AgluonLeg1PlusLimitAware(const BornPhSpPt& born, const RadiationVariables& rad, double muR2) const
    {
        if (rad.y == 1.0)
            return AgluonLeg1PlusCollinearLimit(born, rad, muR2);
        else
        {
            RealPhSpPt real = m_realPhaseSpace->reconstruct(born, rad);
            return AgluonLeg1Plus(real, muR2);
        }
    }

    double BBarIntegrator::AgluonLeg2MinusLimitAware(const BornPhSpPt& born, const RadiationVariables& rad, double muR2) const
    {
        if (rad.y == -1.0)
            return AgluonLeg2MinusCollinearLimit(born, rad, muR2);
        else
        {
            RealPhSpPt real = m_realPhaseSpace->reconstruct(born, rad);
            return AgluonLeg2Minus(real, muR2);
        }
    }

    double BBarIntegrator::RHatgluonLeg1Plus(const BornPhSpPt& born, const std::array<double, 3>& unitCube, double muF2, double muR2) const
    {
        const RadiationVariables rad = sampleUniformRad(born, unitCube);
        const RadiationVariables radColl = sampleUniformRad(born, {unitCube[0], 1.0, unitCube[2]});

        const double lum         = luminosity(born, rad,     21, born.channel.id2, muF2);
        const double lumColl     = luminosity(born, radColl, 21, born.channel.id2, muF2);

        const double A_y = lum     * AgluonLeg1PlusLimitAware(born, rad,     muR2);
        const double A_1 = lumColl * AgluonLeg1PlusLimitAware(born, radColl, muR2);

        return (A_y - A_1) / (1.0 - rad.y);
    }

    double BBarIntegrator::RHatgluonLeg2Minus(const BornPhSpPt& born, const std::array<double, 3>& unitCube, double muF2, double muR2) const
    {
        const RadiationVariables rad = sampleUniformRad(born, unitCube);
        const RadiationVariables radColl = sampleUniformRad(born, {unitCube[0], 0.0, unitCube[2]});

        const double lum         = luminosity(born, rad,     born.channel.id1, 21, muF2);
        const double lumColl     = luminosity(born, radColl, born.channel.id1, 21, muF2);

        const double A_y = lum     * AgluonLeg2MinusLimitAware(born, rad,     muR2);
        const double A_m = lumColl * AgluonLeg2MinusLimitAware(born, radColl, muR2);

        return (A_y - A_m) / (1.0 + rad.y);
    }

    double BBarIntegrator::bornPlusVirtualContribution(const BornPhSpPt& born, double muF2, double muR2) const
    {
        const double f  = m_config.PDF->xfxQ2(born.channel.id1, born.x1Bar, muF2) / born.x1Bar;
        const double fb = m_config.PDF->xfxQ2(born.channel.id2, born.x2Bar, muF2) / born.x2Bar;

        // const double ampBorn = m_process.bornContribution(born);
        // const double ampVirt = m_process.virtualContribution(born, muR2);

        // const double bornContr = born.jacobianOld * f * fb * ampBorn;
        // const double virtContr = born.jacobianOld * f * fb * ampVirt;

        // return bornContr + virtContr;
        return 0.0;
    }

    double BBarIntegrator::realMinusCounterTermContribution(
        const BornPhSpPt& born,
        const std::array<double, 3>& unitCube,
        double muF2,
        double muR2
    ) const 
    {
        const double sigma =
            RHatqqbarPlus(born, unitCube, muF2, muR2) 
            + RHatqqbarMinus(born, unitCube, muF2, muR2)
            + RHatgluonLeg1Plus(born, unitCube, muF2, muR2)
            + RHatgluonLeg2Minus(born, unitCube, muF2, muR2)
            + realEndpointRemnantQQbar(born, unitCube, muF2, muR2)
        ;

        const RadiationVariables rad = sampleUniformRad(born, unitCube);
        const RealPhSpPt real = m_realPhaseSpace->reconstruct(born, rad);

        const double bornPhSpFactor = 1.0 / (64.0 * PI * PI * m_config.S);

        return bornPhSpFactor
            * born.jacobianOld
            * unitCubeJacobian(rad)
            * real.radJacobian
            * sigma;
    }

    double BBarIntegrator::collinearRemnantGluonLeg1(const BornPhSpPt& born, double unitX, double muF2) const
    {
        const double tiny = 0; //1.0e-6; // par_isrtinycsi

        const double z =
            1.0
        - (1.0 - born.x1Bar) * tiny
        - (1.0 - born.x1Bar) * (1.0 - tiny) * unitX;

        const double jac = (1.0 - born.x1Bar) / z;
        
        const double lumGluonLeg1_z = m_config.PDF->xfxQ2(21, born.x1Bar / z, muF2) / (born.x1Bar / z)
            * m_config.PDF->xfxQ2(born.channel.id2, born.x2Bar, muF2) / born.x2Bar;

        const double bracket = (z * z + (1.0 - z) * (1.0 - z)) 
            * (std::log(born.sHat / z / muF2) + 2.0 * std::log(1.0 - z))
            + 2.0 * z * (1.0 - z);

        return m_config.T_F * jac * (lumGluonLeg1_z * bracket);
    }

    double BBarIntegrator::collinearRemnantGluonLeg2(const BornPhSpPt& born, double unitX, double muF2) const
    {
        const double tiny = 0; //1.0e-6; // par_isrtinycsi

        const double z =
            1.0
        - (1.0 - born.x2Bar) * tiny
        - (1.0 - born.x2Bar) * (1.0 - tiny) * unitX;

        const double jac = (1.0 - born.x2Bar) / z;
        
        const double lumGluonLeg2_z = m_config.PDF->xfxQ2(born.channel.id1, born.x1Bar, muF2) / born.x1Bar
            * m_config.PDF->xfxQ2(21, born.x2Bar / z, muF2) / (born.x2Bar / z);

        const double bracket = (z * z + (1.0 - z) * (1.0 - z)) 
            * (std::log(born.sHat / z / muF2) + 2.0 * std::log(1.0 - z))
            + 2.0 * z * (1.0 - z);

        return m_config.T_F * jac * (lumGluonLeg2_z * bracket);
    }

    // double BBarIntegrator::collinearRemnantPlusQQ(const BornPhSpPt& born, double unitX, double muF2) const
    // {
    //     const double z = born.x1Bar + (1.0 - born.x1Bar) * unitX;
    //     const double jac = (1.0 - born.x1Bar) / z;
        
    //     const double lumQQ_z = m_config.PDF->xfxQ2(born.channel.id1, born.x1Bar / z, muF2) / (born.x1Bar / z)
    //         * m_config.PDF->xfxQ2(born.channel.id2, born.x2Bar, muF2) / born.x2Bar;
    //     const double lumQQ_1 = m_config.PDF->xfxQ2(born.channel.id1, born.x1Bar, muF2) / born.x1Bar
    //         * m_config.PDF->xfxQ2(born.channel.id2, born.x2Bar, muF2) / born.x2Bar;

    //     const double distr1 = ((1.0 + z * z) * std::log(born.sHat / z / muF2) * lumQQ_z
    //         - 2.0 * std::log(born.sHat / muF2) * lumQQ_1) / (1.0 - z);
    //     const double distr2 = 2.0 * std::log(1.0 - z) * ((1.0 + z * z) * lumQQ_z - 2.0 * lumQQ_1) / (1.0 - z);

    //     return m_config.C_F * jac * (distr1 + distr2 + lumQQ_z * (1.0 - z));
    // }

    double BBarIntegrator::collinearRemnantPlusQQ(
        const BornPhSpPt& born,
        double unitX,
        double muF2
    ) const
    {
        const double tiny = 0; //1.0e-6; // par_isrtinycsi

        const double z = 1.0 - (1.0 - born.x1Bar) * tiny - (1.0 - born.x1Bar) * (1.0 - tiny) * unitX;
        const double dz = (1.0 - born.x1Bar) * (1.0 - tiny);

        const double jacOverZ = dz / z;
        const double jacAtOne = dz;

        const double lumQQ_z = m_config.PDF->xfxQ2(born.channel.id1, born.x1Bar / z, muF2) / (born.x1Bar / z)
            * m_config.PDF->xfxQ2(born.channel.id2, born.x2Bar, muF2) / born.x2Bar;

        const double lumQQ_1 = m_config.PDF->xfxQ2(born.channel.id1, born.x1Bar, muF2) / born.x1Bar
            * m_config.PDF->xfxQ2(born.channel.id2, born.x2Bar, muF2) / born.x2Bar;

        const double logReal = std::log(born.sHat / z / muF2) + 2.0 * std::log(1.0 - z);
        const double logSub  = std::log(born.sHat / muF2)     + 2.0 * std::log(1.0 - z);

        const double realScaled =
            jacOverZ * lumQQ_z *
            ((1.0 + z * z) * logReal / (1.0 - z) + (1.0 - z));

        const double subtractionAtOne =
            jacAtOne * lumQQ_1 *
            (2.0 * logSub / (1.0 - z));

        const double csicut = 1.0; // match POWHEG par_csicut

        const double rm1 = std::log((1.0 - born.x1Bar) / csicut) * std::log(born.sHat / muF2)
            + std::pow(std::log(1.0 - born.x1Bar), 2)
            - std::pow(std::log(csicut), 2);

        const double deltaRemnant = 2.0 * lumQQ_1 * rm1;

        return m_config.C_F * (realScaled - subtractionAtOne 
            + deltaRemnant
        );
    }

    // double BBarIntegrator::collinearRemnantMinusQQ(const BornPhSpPt& born, double unitX, double muF2) const
    // {
    //     const double z = born.x2Bar + (1.0 - born.x2Bar) * unitX;
    //     const double jac = (1.0 - born.x2Bar) / z;
        
    //     const double lumQQ_z = m_config.PDF->xfxQ2(born.channel.id1, born.x1Bar, muF2) / born.x1Bar
    //         * m_config.PDF->xfxQ2(born.channel.id2, born.x2Bar / z, muF2) / (born.x2Bar / z);
    //     const double lumQQ_1 = m_config.PDF->xfxQ2(born.channel.id1, born.x1Bar, muF2) / born.x1Bar
    //         * m_config.PDF->xfxQ2(born.channel.id2, born.x2Bar, muF2) / born.x2Bar;

    //     const double distr1 = ((1.0 + z * z) * std::log(born.sHat / z / muF2) * lumQQ_z
    //         - 2.0 * std::log(born.sHat / muF2) * lumQQ_1) / (1.0 - z);
    //     const double distr2 = 2.0 * std::log(1.0 - z) * ((1.0 + z * z) * lumQQ_z - 2.0 * lumQQ_1) / (1.0 - z);

    //     return m_config.C_F * jac * (distr1 + distr2 + lumQQ_z * (1.0 - z));
    // }

    double BBarIntegrator::collinearRemnantMinusQQ(
        const BornPhSpPt& born,
        double unitX,
        double muF2
    ) const
    {
        const double tiny = 0; //1.0e-6; // par_isrtinycsi

        const double z =
            1.0
        - (1.0 - born.x2Bar) * tiny
        - (1.0 - born.x2Bar) * (1.0 - tiny) * unitX;

        const double dz = (1.0 - born.x2Bar) * (1.0 - tiny);

        const double jacOverZ = dz / z;
        const double jacAtOne = dz;

        const double lumQQ_z =
            m_config.PDF->xfxQ2(born.channel.id1, born.x1Bar, muF2) / born.x1Bar
        * m_config.PDF->xfxQ2(born.channel.id2, born.x2Bar / z, muF2) / (born.x2Bar / z);

        const double lumQQ_1 =
            m_config.PDF->xfxQ2(born.channel.id1, born.x1Bar, muF2) / born.x1Bar
        * m_config.PDF->xfxQ2(born.channel.id2, born.x2Bar, muF2) / born.x2Bar;

        const double logReal = std::log(born.sHat / z / muF2) + 2.0 * std::log(1.0 - z);
        const double logSub  = std::log(born.sHat / muF2)     + 2.0 * std::log(1.0 - z);

        const double realScaled =
            jacOverZ * lumQQ_z *
            ((1.0 + z * z) * logReal / (1.0 - z) + (1.0 - z));

        const double subtractionAtOne =
            jacAtOne * lumQQ_1 *
            (2.0 * logSub / (1.0 - z));

        const double csicut = 1.0; // match POWHEG par_csicut

        const double rm2 = std::log((1.0 - born.x2Bar) / csicut) * std::log(born.sHat / muF2)
            + std::pow(std::log(1.0 - born.x2Bar), 2)
            - std::pow(std::log(csicut), 2);

        const double deltaRemnant = 2.0 * lumQQ_1 * rm2;

        return m_config.C_F * (realScaled - subtractionAtOne
             + deltaRemnant
            );
    }

    double BBarIntegrator::realEndpointRemnantQQbar(
        const BornPhSpPt& born,
        const std::array<double, 3>& unitCube,
        double muF2,
        double muR2
    ) const
    {
        const double csiCut = 1.0; // Use POWHEG par_csicut here.

        const RadiationVariables rad = sampleUniformRad(born, unitCube);
        const RadiationVariables radSoft = sampleUniformRad(born, {0.0, unitCube[1], unitCube[2]});
        const RadiationVariables radSoftCollPlus = sampleUniformRad(born, {0.0, 1.0, unitCube[2]});
        const RadiationVariables radSoftCollMinus = sampleUniformRad(born, {0.0, 0.0, unitCube[2]});

        const double xi = rad.xi;
        const double y  = rad.y;

        const double xiTilde = xi / rad.xiMax;

        const double lumSoft = luminosity(born, radSoft, born.channel.id1, born.channel.id2, muF2);
        const double lumSoftCollPlus = luminosity(born, radSoftCollPlus, born.channel.id1, born.channel.id2, muF2);
        const double lumSoftCollMinus = luminosity(born, radSoftCollMinus, born.channel.id1, born.channel.id2, muF2);

        const double A_0y_plus = lumSoft * AqqbarPlusSoftLimit(born, radSoft, muR2);
        const double A_0y_minus = lumSoft * AqqbarMinusSoftLimit(born, radSoft, muR2);
        const double A_01_plus = lumSoftCollPlus * AqqbarPlusSoftCollinearLimit(born, muR2);
        const double A_0m_minus = lumSoftCollMinus * AqqbarMinusSoftCollinearLimit(born, muR2);

        const double xi2 = xi * xi;

        const double soft = A_0y_plus  / (xi2 * (1.0 - y)) + A_0y_minus / (xi2 * (1.0 + y));

        const double softCollPlus = 0.5 * A_01_plus / (xi2 * (1.0 - y));
        const double softCollMinus = 0.5 * A_0m_minus / (xi2 * (1.0 + y));

        const double logSoft = std::log(rad.xiMax / csiCut);
        const double logSoftCollPlus = std::log((1.0 - born.x1Bar) / csiCut);
        const double logSoftCollMinus = std::log((1.0 - born.x2Bar) / csiCut);

        const double softOverRealJac = (1.0 - xi) * (1.0 - xi);

        return softOverRealJac * xiTilde * (
            soft          * logSoft
            - softCollPlus  * logSoftCollPlus
            - softCollMinus * logSoftCollMinus
        );
    }

    double BBarIntegrator::collinearRemnantContribution(const BornPhSpPt& born, const std::array<double, 4>& unitX, double muF2, double muR2) const
    {
        const double prefactor = 1.0 / (64.0 * PI * PI * m_config.S * born.sHat);

        // TODO: Maybe go back to using all for unitX values
        const double plusQQ = collinearRemnantPlusQQ(born, unitX[0], muF2);
        const double minusQQ = collinearRemnantMinusQQ(born, unitX[0], muF2);
        const double gluonLeg1 = collinearRemnantGluonLeg1(born, unitX[0], muF2);
        const double gluonLeg2 = collinearRemnantGluonLeg2(born, unitX[0], muF2);

        return prefactor * born.jacobianOld * m_config.alphaS(muR2) / 2.0 / PI * (plusQQ + minusQQ + gluonLeg1 + gluonLeg2)
            * m_process.bornAmp2(born);
    }

    double BBarIntegrator::bTilde(
        const BornPhSpPt& born,
        const std::array<double, 3>& unitCube
    ) const
    {
        const double muF2 = born.sHat;
        const double muR2 = born.sHat;

        const BBarIntegrationPoint point = generateIntegrationPoint(born, muF2, muR2, unitCube);

        double dSigma = m_bornVirtual.dSigmaBorn(point);

        if (m_config.BTILDE_BORNONLY)
        return dSigma;
        
        dSigma += m_bornVirtual.dSigmaVirtual(point, muR2);

        const double dSigmaBornVirt = dSigma;

        for (const auto& realChannel : m_process.realChannels(born.channel))
        {
            for (const auto& singularRegion : findSingularRegions(born.channel, realChannel))
            {
                if (const auto& collinearChannel = remnantChannelFromRegion(born.channel, realChannel, singularRegion); 
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

        // const double dSigmaColl = dSigma - dSigmaBornVirt;
        // const double absTot = std::abs(dSigmaBornVirt) + std::abs(dSigmaColl);
        // p1 += std::abs(dSigmaBornVirt) / absTot;
        // p2 += std::abs(real) / absTot;
        // p3 += std::abs(dSigmaColl) / absTot;
        // n++;

        return dSigma;
    }

    double BBarIntegrator::bTildeOld(BornPhSpPt& born, const std::array<double, 3>& unitCube, const std::array<double, 4>& unitX) const
    {
        const double muF2 = born.sHat;
        const double muR2 = born.sHat;

        // TODO: Move this block into a different function??
        born.f1 = m_config.PDF->xfxQ2(born.channel.id1, born.x1Bar, muF2) / born.x1Bar;
        born.f2 = m_config.PDF->xfxQ2(born.channel.id2, born.x2Bar, muF2) / born.x2Bar;
        born.amp2 = m_process.bornAmp2(born);

        // const double bornPlusVirtual = bornPlusVirtualContribution(born, muF2, muR2);
        // const double realMinusCounterterm = realMinusCounterTermContribution(born, unitCube, muF2, muR2);
        // const double collinearRemnants = collinearRemnantContribution(born, unitX, muF2, muR2);

        const double absTot = 0.0
            // + std::abs(dSigmaBorn + dSigmaVirtual) 
            // + std::abs(realMinusCounterterm)
            // + std::abs(collinearRemnants)
        ;

        // p1 += std::abs(dSigmaBorn + dSigmaVirtual) / absTot;
        // p2 += std::abs(realMinusCounterterm) / absTot;
        // p3 += std::abs(collinearRemnants) / absTot;
        n++;

        return 0.0
            // + dSigmaBorn 
            // + dSigmaVirtual 
            // + realMinusCounterterm
            // + collinearRemnants
        ;
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
                double rands[3] = { rand(), rand(), rand() }; 

                BornPhSpPt born = m_bornPhaseSpace->samplePoint(rands);
                computeWeightAndSampleChannel(born);

                if (born.weight > maxWeight)
                    maxWeight = born.weight;
            }

            m_maxWeight = SECURITY_FACTOR * maxWeight;

            Log::info << "Done determining Born Veto weight." << std::endl << std::endl;
        }
    }

    double BBarIntegrator::getTotalCrossSection() const 
    {
        static int i = 0;
        if (++i == 1)
        {
            Log::info << "Born average relative contr. " << p1 / n << std::endl;
            Log::info << "Real average relative contr. " << p2 / n << std::endl;
            Log::info << "Coll average relative contr. " << p3 / n << std::endl;
            Log::info << "Number of negative events: " << (m_config.N_ACCEPTED_EVENTS - m_sumSigns) / 2 << std::endl;
            Log::info << "Veto weight: " << m_maxWeight << std::endl;
        }

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
        const double muF2,
        const double muR2,
        const std::array<double, 3>& unitCube
    ) const
    {
        BBarIntegrationPoint point = BBarIntegrationPoint{ born };

        point.u1 = unitCube[0];
        point.u2 = unitCube[1];
        point.u3 = unitCube[2];

        point.y = -1.0 + 2.0 * point.u2;
        point.xi = m_realPhaseSpace->xiMax(born, point.y) * point.u1;
        point.phi = 2.0 * PI * point.u3;

        point.unitCubeJacobian = 4.0 * PI * point.xi / point.u1;        

        point.zLeg1 = born.x1Bar + (1.0 - born.x1Bar) * point.u1;
        point.zLeg2 = born.x2Bar + (1.0 - born.x2Bar) * point.u1;

        point.f1Born = m_config.PDF->xfxQ2(born.channel.id1, born.x1Bar, muF2) / born.x1Bar; 
        point.f2Born = m_config.PDF->xfxQ2(born.channel.id2, born.x2Bar, muF2) / born.x2Bar; 

        #warning probably these real luminosities are only for collinear remnants, not for the real part
        point.f1RealQ = m_config.PDF->xfxQ2(born.channel.id1, born.x1Bar / point.zLeg1, muF2) / (born.x1Bar / point.zLeg1); 
        point.f2RealQ = m_config.PDF->xfxQ2(born.channel.id2, born.x2Bar / point.zLeg2, muF2) / (born.x2Bar / point.zLeg2); 

        point.f1RealG = m_config.PDF->xfxQ2(21, born.x1Bar / point.zLeg1, muF2) / (born.x1Bar / point.zLeg1); 
        point.f2RealG = m_config.PDF->xfxQ2(21, born.x2Bar / point.zLeg2, muF2) / (born.x2Bar / point.zLeg2); 
        
        point.alphaS = m_config.alphaS(muR2);

        point.amp2Born = m_process.bornAmp2(born);

        return point;
    }

} // namespace powheg_dy
