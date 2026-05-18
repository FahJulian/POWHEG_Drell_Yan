#pragma once

#include "config.h"
#include "powheg/process.h"

namespace powheg
{
    class DrellYanProcess : public Process<DrellYanConfig>
    {
        void initConfig(ConfigParser& parser) const override;

        double bornAmp2(const BornPhSpPt& born, const BornChannel& bornChannel) const override;
        double virtAmp2(const BornPhSpPt& born, const BornChannel& bornChannel, double amp2Born, double muR2) const override;
        double realAmp2(const RealPhSpPt& real, const RealChannel& realChannel, double alphaS) const override;

        std::vector<BornChannel> bornChannels() const override;
        std::vector<RealChannel> realChannels(const BornChannel& bornChannel) const override;

        BornPhSpPt sampleBorn() const override;

        double bornAmp2(const int flavour, const FourVector& pQ, const FourVector& pQbar, 
            const FourVector& pLMinus, const FourVector& pLPlus) const;
        double realAmp2_crossed(const int flavour, const FourVector& qSlot, const FourVector& qbarSlot, const FourVector& gluonSlot, 
            const FourVector& pLMinus, const FourVector& pLPlus, double alphaS, double colorFactor) const;

        double realAmp2qqbar(const RealPhSpPt& real, const RealChannel& realChannel, double alphaS) const;
        double realAmp2gluonLeg1(const RealPhSpPt& real, const RealChannel& realChannel, double alphaS) const;
        double realAmp2gluonLeg2(const RealPhSpPt& real, const RealChannel& realChannel, double alphaS) const;

        std::pair<double, double> sampleMSq() const;
        std::pair<double, double> sampleYBoson(const double sHat) const;
        std::pair<double, double> sampleCosTh() const;
        std::pair<double, double> samplePhi() const;
    };

} // namespace powheg
