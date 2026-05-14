#pragma once

#include "powheg_dy/process.h"

namespace powheg_dy
{
    class DrellYanProcess : public Process
    {
    public:
        enum RealChannelOld : int
        {
            QQBAR = 0,
            GLUON_LEG1 = 1,
            GLUON_LEG2 = 2,
        };

        double bornAmp2(const BornPhSpPt& born) const override;
        double virtAmp2(const BornPhSpPt& born, const double amp2Born, const double muR2) const override;
        double realAmp2(const RealPhSpPt& real, const RealChannel& channel, const double alphaS) const override;

        std::vector<BornChannel> bornChannels() const override;
        std::vector<RealChannel> realChannels(const BornChannel& bornChannel) const override;

    private:        
        double bornAmp2(int flavour, const FourVector& pQ, const FourVector& pQbar, const FourVector& pLMinus, 
                const FourVector& pLPlus) const;

        double realAmp2_crossed( int flavour, const FourVector& qSlot, const FourVector& qbarSlot, const FourVector& gluonSlot, 
                const FourVector& pLMinus, const FourVector& pLPlus, double alphaS, double colorFactor) const;

        double realAmp2qqbar(const RealPhSpPt& real, const double alphaS) const;
        double realAmp2gluonLeg1(const RealPhSpPt& real, const double alphaS) const;
        double realAmp2gluonLeg2(const RealPhSpPt& real, const double alphaS) const;
    };

} // namespace powheg_dy
