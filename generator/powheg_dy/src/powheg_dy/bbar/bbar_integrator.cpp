#include "bbar_integrator.h"

#include "powheg_dy/process.h"
#include "powheg_dy/alpha_s.h"
#include "powheg_dy/math/rand.h"

namespace powheg_dy
{
namespace 
{
    static constexpr double GEV_SQ_TO_PB = 0.389379338e9;
    static constexpr double SECURITY_FACTOR = 1.1;
    static constexpr int MAX_TRIALS = 10000;

    double sPlus(double y)
    {
        return 0.5 * (1.0 + y);
    }

    double sMinus(double y) 
    {
        return 0.5 * (1.0 - y);
    }

    enum RealChannel : int
    {
        QQBAR = 0,
        GLUON_LEG1 = 1,
        GLUON_LEG2 = 2,
    };

} // anonymous namespace 
    
    bool BBarIntegrator::computeWeightAndSampleChannel(BornPhSpPt& born) const
    {
        const std::array<double, 3> unitCube = { rand(), rand(), rand() };      // For the radiation variables

        std::vector<BornChannel> channels = m_process.bornChannels();
        
        std::vector<std::pair<BornPhSpPt, double>> weights = {};
        weights.reserve(channels.size());
            
        double totalWeight = 0.0;

        for (const auto& channel : channels)
        {
            BornPhSpPt bornCopy = born;
            bornCopy.channel = channel;
            m_bornPhaseSpace->reconstructMomenta(bornCopy);

            const double weight = bTilde(bornCopy, unitCube);

            totalWeight += weight;
            weights.push_back({ bornCopy, weight });
        }

        // assert(totalWeight > 0.0);
        if (totalWeight <= 0.0)
        {
            // Log::warn("Negative Weight!");
            return false;
        }

        born.weight = totalWeight;

        // Sample the parton channel by their relative contribution to dSigma
        double u = rand(0.0, totalWeight);
        for (auto& [bornCopy, weight] : weights)
        {
            if (u < weight)
            {   
                bornCopy.weight = totalWeight;
                born = bornCopy;

                break;
            }

            u -= weight;
        }

        return true;
    }

    BornPhSpPt BBarIntegrator::sampleAccordingtoBTilde()
    {
        for (int trials = 1; trials <= MAX_TRIALS; trials++)
        {
            double rands[3] = { rand(), rand(), rand() };
            BornPhSpPt born = m_bornPhaseSpace->samplePoint(rands);
            
            if (!computeWeightAndSampleChannel(born))
                continue;

            assert(born.weight <= m_maxWeight, "Born weight " << born.weight << " exceeds max weight " << m_maxWeight);

            double u = rand();
            if (u < born.weight / m_maxWeight)
            {
                m_nEventTrials += trials;
                return born;
            }
        }

        assert(false, "Bbar sampling stuck in loop.");
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

        // if (x1 <= 0.0 || x1 >= 1.0 || x2 <= 0.0 || x2 >= 1.0)
        // {
        //     Log::warn << "Unphysical values of x: x1 = " << x1 << ", x2 = " << x2 << std::endl;
        //     return 0.0;
        // }
        
        const double fLeg1 = m_config.PDF->xfxQ2(id1, x1, muF2) / x1; 
        const double fLeg2 = m_config.PDF->xfxQ2(id2, x2, muF2) / x2; 

        return fLeg1 * fLeg2;
    }

    double BBarIntegrator::AqqbarPlus(const RealPhSpPt& real, double muF2, double muR2) const
    {
        return sPlus(real.rad.y) * (1.0 - real.rad.y) * real.rad.xi * real.rad.xi 
            * m_process.realAmp2(real, muR2, RealChannel::QQBAR, false) / real.sHatReal;
    }

    double BBarIntegrator::AqqbarMinus(const RealPhSpPt& real, double muF2, double muR2) const
    {
        return sMinus(real.rad.y) * (1.0 + real.rad.y) * real.rad.xi * real.rad.xi 
            * m_process.realAmp2(real, muR2, RealChannel::QQBAR, false) / real.sHatReal;
    }

    double BBarIntegrator::AqqbarPlusLimitAware(const BornPhSpPt& born, const RadiationVariables& rad, double muF2, double muR2) const
    {
        const bool soft = (rad.xi == 0.0);
        const bool coll = (rad.y  == 1.0);

        if (soft && coll)
            return AqqbarPlusSoftCollinearLimit(born, rad, muF2, muR2);
        else if (soft)
            return AqqbarPlusSoftLimit(born, rad, muF2, muR2);
        else if (coll)
            return AqqbarPlusCollinearLimit(born, rad, muF2, muR2);
        else
        {
            RealPhSpPt real = m_realPhaseSpace->reconstruct(born, rad);
            return AqqbarPlus(real, muF2, muR2);
        }
    }

    double BBarIntegrator::AqqbarMinusLimitAware(const BornPhSpPt& born, const RadiationVariables& rad, double muF2, double muR2) const
    {
        const bool soft = (rad.xi ==  0.0);
        const bool coll = (rad.y  == -1.0);

        if (soft && coll)
            return AqqbarMinusSoftCollinearLimit(born, rad, muF2, muR2);
        else if (soft)
            return AqqbarMinusSoftLimit(born, rad, muF2, muR2);
        else if (coll)
            return AqqbarMinusCollinearLimit(born, rad, muF2, muR2);
        else
        {
            RealPhSpPt real = m_realPhaseSpace->reconstruct(born, rad);
            return AqqbarMinus(real, muF2, muR2);
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

        const double A_xy = lum         * AqqbarPlusLimitAware(born, rad,         muF2, muR2);
        const double A_0y = lumSoft     * AqqbarPlusLimitAware(born, radSoft,     muF2, muR2);
        const double A_x1 = lumColl     * AqqbarPlusLimitAware(born, radColl,     muF2, muR2);
        const double A_01 = lumSoftColl * AqqbarPlusLimitAware(born, radSoftColl, muF2, muR2);

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

        const double A_xy = lum         * AqqbarMinusLimitAware(born, rad,         muF2, muR2);
        const double A_0y = lumSoft     * AqqbarMinusLimitAware(born, radSoft,     muF2, muR2);
        const double A_xm = lumColl     * AqqbarMinusLimitAware(born, radColl,     muF2, muR2);
        const double A_0m = lumSoftColl * AqqbarMinusLimitAware(born, radSoftColl, muF2, muR2);

        return (A_xy - A_0y - A_xm + A_0m) / (rad.xi * rad.xi * (1.0 + rad.y));
    }

    double BBarIntegrator::AqqbarPlusCollinearLimit(const BornPhSpPt& born, const RadiationVariables& rad, double muF2, double muR2) const
    {
        const double z = 1.0 - rad.xi;
        const double prefactor = 16.0 * PI * m_config.C_F / born.sHat / born.sHat;

        return prefactor
            * alphaS(m_config, muR2)
            * m_process.bornAmp2(born)
            * z * (1 + z * z);
    }

    double BBarIntegrator::AqqbarMinusCollinearLimit(const BornPhSpPt& born, const RadiationVariables& rad, double muF2, double muR2) const
    {
        const double z = 1.0 - rad.xi;
        const double prefactor = 16.0 * PI * m_config.C_F / born.sHat / born.sHat;

        return prefactor
            * alphaS(m_config, muR2)
            * m_process.bornAmp2(born)
            * z * (1 + z * z);
    }

    double BBarIntegrator::AqqbarPlusSoftLimit(const BornPhSpPt& born, const RadiationVariables& rad, double muF2, double muR2) const
    {
        const double angular = (1.0 - rad.y) * sPlus(rad.y) / (1.0 - rad.y * rad.y);
        const double prefactor = 64.0 * PI * m_config.C_F / born.sHat / born.sHat;

        return prefactor
            * alphaS(m_config, muR2)
            * m_process.bornAmp2(born)
            * angular;
    }

    double BBarIntegrator::AqqbarMinusSoftLimit(const BornPhSpPt& born, const RadiationVariables& rad, double muF2, double muR2) const
    {
        const double angular = (1.0 + rad.y) * sMinus(rad.y) / (1.0 - rad.y * rad.y);
        const double prefactor = 64.0 * PI * m_config.C_F / born.sHat / born.sHat;

        return prefactor
            * alphaS(m_config, muR2)
            * m_process.bornAmp2(born)
            * angular;
    }

    double BBarIntegrator::AqqbarPlusSoftCollinearLimit(const BornPhSpPt& born, const RadiationVariables& rad, double muF2, double muR2) const
    {
        const double prefactor = 32.0 * PI * m_config.C_F / born.sHat / born.sHat;
        return prefactor * alphaS(m_config, muR2) * m_process.bornAmp2(born);
    }

    double BBarIntegrator::AqqbarMinusSoftCollinearLimit(const BornPhSpPt& born, const RadiationVariables& rad, double muF2, double muR2) const
    {
        const double prefactor = 32.0 * PI * m_config.C_F / born.sHat / born.sHat;
        return prefactor * alphaS(m_config, muR2) * m_process.bornAmp2(born);
    }

    double BBarIntegrator::AgluonLeg1Plus(const RealPhSpPt& real, double muF2, double muR2) const
    {
        return (1.0 - real.rad.y) * m_process.realAmp2(real, muR2, RealChannel::GLUON_LEG1, false) / real.sHatReal;
    }

    double BBarIntegrator::AgluonLeg2Minus(const RealPhSpPt& real, double muF2, double muR2) const
    {
        return (1.0 + real.rad.y) * m_process.realAmp2(real, muR2, RealChannel::GLUON_LEG2, false) / real.sHatReal;
    }

    double BBarIntegrator::AgluonLeg1PlusCollinearLimit(const BornPhSpPt& born, const RadiationVariables& rad, double muF2, double muR2) const
    {
        const double z = 1.0 - rad.xi;
        const double prefactor = 16.0 * PI * m_config.T_F / born.sHat / born.sHat;

        return m_process.bornAmp2(born)
            * prefactor
            * alphaS(m_config, muR2)
            * z
            * (z * z + (1.0 - z) * (1.0 - z))
            / (1.0 - z);
    }

    double BBarIntegrator::AgluonLeg2MinusCollinearLimit(const BornPhSpPt& born, const RadiationVariables& rad, double muF2, double muR2) const
    {
        const double z = 1.0 - rad.xi;
        const double prefactor = 16.0 * PI * m_config.T_F / born.sHat / born.sHat;

        return m_process.bornAmp2(born)
            * prefactor
            * alphaS(m_config, muR2)
            * z
            * (z * z + (1.0 - z) * (1.0 - z))
            / (1.0 - z);
    }

    double BBarIntegrator::AgluonLeg1PlusLimitAware(const BornPhSpPt& born, const RadiationVariables& rad, double muF2, double muR2) const
    {
        if (rad.y == 1.0)
            return AgluonLeg1PlusCollinearLimit(born, rad, muF2, muR2);
        else
        {
            RealPhSpPt real = m_realPhaseSpace->reconstruct(born, rad);
            return AgluonLeg1Plus(real, muF2, muR2);
        }
    }

    double BBarIntegrator::AgluonLeg2MinusLimitAware(const BornPhSpPt& born, const RadiationVariables& rad, double muF2, double muR2) const
    {
        if (rad.y == -1.0)
            return AgluonLeg2MinusCollinearLimit(born, rad, muF2, muR2);
        else
        {
            RealPhSpPt real = m_realPhaseSpace->reconstruct(born, rad);
            return AgluonLeg2Minus(real, muF2, muR2);
        }
    }

    double BBarIntegrator::RHatgluonLeg1Plus(const BornPhSpPt& born, const std::array<double, 3>& unitCube, double muF2, double muR2) const
    {
        const RadiationVariables rad = sampleUniformRad(born, unitCube);
        const RadiationVariables radColl = sampleUniformRad(born, {unitCube[0], 1.0, unitCube[2]});

        const double lum         = luminosity(born, rad,     21, born.channel.id2, muF2);
        const double lumColl     = luminosity(born, radColl, 21, born.channel.id2, muF2);

        const double A_y = lum     * AgluonLeg1PlusLimitAware(born, rad,     muF2, muR2);
        const double A_1 = lumColl * AgluonLeg1PlusLimitAware(born, radColl, muF2, muR2);

        return (A_y - A_1) / (1.0 - rad.y);
    }

    double BBarIntegrator::RHatgluonLeg2Minus(const BornPhSpPt& born, const std::array<double, 3>& unitCube, double muF2, double muR2) const
    {
        const RadiationVariables rad = sampleUniformRad(born, unitCube);
        const RadiationVariables radColl = sampleUniformRad(born, {unitCube[0], 0.0, unitCube[2]});

        const double lum         = luminosity(born, rad,     born.channel.id1, 21, muF2);
        const double lumColl     = luminosity(born, radColl, born.channel.id1, 21, muF2);

        const double A_y = lum     * AgluonLeg2MinusLimitAware(born, rad,     muF2, muR2);
        const double A_m = lumColl * AgluonLeg2MinusLimitAware(born, radColl, muF2, muR2);

        return (A_y - A_m) / (1.0 + rad.y);
    }

    BBarIntegrator::Counterterms BBarIntegrator::counterterms(const BornPhSpPt& born, const RadiationVariables& rad) const
    {
        return { };
    }

    double BBarIntegrator::bornPlusVirtualContribution(const BornPhSpPt& born) const
    {
        const double muF2 = born.sHat;
        const double muR2 = born.sHat;

        const double f  = m_config.PDF->xfxQ2(born.channel.id1, born.x1Bar, muF2) / born.x1Bar;
        const double fb = m_config.PDF->xfxQ2(born.channel.id2, born.x2Bar, muF2) / born.x2Bar;

        const double ampBorn = m_process.bornContribution(born);
        const double ampVirt = m_process.virtualContribution(born, muR2);

        const double bornContr = born.jacobian * f * fb * ampBorn;
        const double virtContr = born.jacobian * f * fb * ampVirt;

        return bornContr + virtContr;
    }

    double BBarIntegrator::realMinusCounterTermContribution(const BornPhSpPt& born, const std::array<double, 3>& unitCube) const 
    {
        const double muF2 = born.sHat;
        const double muR2 = born.sHat;

        // RealPhSpPt real = m_realPhaseSpace->reconstruct(born, rad);
        
        // const double fReal1Q = m_config.PDF->xfxQ2(born.channel.id1, real.x1, muF2) / real.x1; 
        // const double fReal2Q = m_config.PDF->xfxQ2(born.channel.id2, real.x2, muF2) / real.x2; 

        // const double fReal1G = m_config.PDF->xfxQ2(21, real.x1, muF2) / real.x1; 
        // const double fReal2G = m_config.PDF->xfxQ2(21, real.x2, muF2) / real.x2; 

        // const double lumQQbar = fReal1Q * fReal2Q;
        // const double lumGLeg1 = fReal1G * fReal2Q;
        // const double lumGLeg2 = fReal1Q * fReal2G;

        const double sigma = RHatqqbarPlus(born, unitCube, muF2, muR2) 
            + RHatqqbarMinus(born, unitCube, muF2, muR2)
            + RHatgluonLeg1Plus(born, unitCube, muF2, muR2)
            + RHatgluonLeg2Minus(born, unitCube, muF2, muR2);
        
        // TODO: Real is getting sampled too often
        auto rad = sampleUniformRad(born, unitCube);
        auto real = m_realPhaseSpace->reconstruct(born, rad);

        const double bornPhSpFactor = 1.0 / (64.0 * PI * PI * m_config.S);

        return bornPhSpFactor * born.jacobian * unitCubeJacobian(rad) * real.radJacobian * sigma;
    }

    double BBarIntegrator::bTilde(const BornPhSpPt& born, const std::array<double, 3>& unitCube) const
    {
        const double bornPlusVirtual = bornPlusVirtualContribution(born);
        const double realMinusCounterterm = realMinusCounterTermContribution(born, unitCube);

        return 0.0
            + bornPlusVirtual 
            + realMinusCounterterm
        ;
    }

    void BBarIntegrator::determineMaxWeight()
    {
        Log::info("Determining Born Veto weight");

        double max_dSigma = 0.0;
        
        for (size_t i = 0; i < m_config.N_TRIAL_EVENTS; i++)
        {   
            double rands[3] = { rand(), rand(), rand() }; 

            BornPhSpPt born = m_bornPhaseSpace->samplePoint(rands);
            if (!computeWeightAndSampleChannel(born))
            {
                i -= 1;
                continue;
            }

            if (born.weight > max_dSigma)
                max_dSigma = born.weight;
        }

        m_maxWeight = SECURITY_FACTOR * max_dSigma;

        Log::info << "Done determining Born Veto weight." << std::endl << std::endl;
    }

    double BBarIntegrator::getTotalCrossSection() const 
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
    }

    // void test()
    // {
        // auto rad2 = rad;

        // born.channel = { -1, 1, 1 };
        // m_bornPhaseSpace->reconstructMomenta(born);

        // Log::info << "--------------- GLUON LEG1 PLUS COLLINEAR LIMIT --------------" << std::endl;

        // auto [rad2, radJacobian2] = sampleUniformRad(born);
        // for (double n = 1e3; n <= 1e10; n *= 10)
        // {
        //     rad2.y = 1.0 - 1.0 / n;

        //     const double a1 = AgluonLeg1Plus(born, rad2, born.sHat);
        //     const double a2 = AgluonLeg1PlusCollinearLimit(born, rad2.xi, born.sHat);

        //     Log::info << "n = " << n << ": " << a1 / a2 << std::endl;
        // }

        // Log::info << "--------------- GLUON LEG2 MINUS COLLINEAR LIMIT --------------" << std::endl;

        // auto [rad2, radJacobian2] = sampleUniformRad(born);
        // for (double n = 1e3; n <= 1e10; n *= 10)
        // {
        //     rad2.y = -1.0 + 1.0 / n;

        //     const double a1 = AgluonLeg2Minus(born, rad2, born.sHat);
        //     const double a2 = AgluonLeg2MinusCollinearLimit(born, rad2.xi, born.sHat);

        //     Log::info << "n = " << n << ": " << a1 / a2 << std::endl;
        // }

        // Log::info << "--------------- SOFT LIMIT --------------" << std::endl;
        // for (int j = 0; j < 10; j++)
        // {
        //     auto [rad2, radJacobian2] = sampleUniformRad(born);
        //     for (double n = 1e7; n <= 1e10; n *= 10)
        //     {  
        //         rad2.xi = 1.0 / n;
        //         double a1 = AqqbarMinus(born, rad2, born.sHat);
        //         double a2 = AqqbarMinusSoftLimit(born, rad2.y, born.sHat);
        //         Log::info << "n = " << n << ": " << a1 / a2 << std::endl;
        //     }
        // }

        // Log::info << "--------------- COLLINEAR LIMIT --------------" << std::endl;
        // auto [rad2, radJacobian2] = sampleUniformRad(born);
        // for (int j = 0; j < 10; j++)
        // {
        //     double rands[3] = { rand(), rand(), rand() };
        //     auto born2 = m_bornPhaseSpace->samplePoint(rands);
        //     born2.channel = { -1, 1, 1 };
        //     m_bornPhaseSpace->reconstructMomenta(born2);

        //     for (double n = 1e5; n <= 1e10; n *= 10)
        //     {  
        //         rad2.y = 1.0 - 1.0 / n;
        //         double a1 = AqqbarPlus(born2, rad2, born2.sHat);
        //         double a2 = AqqbarPlusCollinearLimit(born2, rad2.xi, born2.sHat);
        //         Log::info << "n = " << n << ": " << a1 / a2 << std::endl;
        //     }
        // }

        // Log::info << "--------------- SOFT COLLINEAR LIMIT --------------" << std::endl;
        // auto [rad2, radJacobian2] = sampleUniformRad(born);
        // for (int j = 0; j < 10; j++)
        // {
        //     double rands[3] = { rand(), rand(), rand() };
        //     auto born2 = m_bornPhaseSpace->samplePoint(rands);
        //     born2.channel = { -1, 1, 1 };
        //     m_bornPhaseSpace->reconstructMomenta(born2);

        //     for (double n = 1e5; n <= 1e10; n *= 10)
        //     {  
        //         rad2.xi = 1.0 / n;
        //         rad2.y  = 1.0 - 1.0 / n;
        //         double a1 = AqqbarPlus(born2, rad2, born2.sHat);
        //         double a2 = AqqbarPlusSoftCollinearLimit(born2, born2.sHat);
        //         Log::info << "n = " << n << ": " << a1 / a2 << std::endl;
        //     }
        // }
    // }

} // namespace powheg_dy
