#pragma once

#include "powheg_dy/process.h"

namespace powheg_dy
{
    class DrellYanProcess : public Process
    {
    public:
        enum RealChannel : int
        {
            QQBAR = 0,
            GLUON_LEG1 = 1,
            GLUON_LEG2 = 2,
        };

        double bornAmp2(const BornPhSpPt& born) const override;

        double realAmp2(const RealPhSpPt real, double muR2, int channel, bool useCMWAlphaS) const override;

        RealOverBornContributions realOverBornContributions(
            const RealPhSpPt& real, double muF2, double muR2, bool useCMWALphaS) const override;

        double bornContribution(const BornPhSpPt& born) const override;
        double virtualContribution(const BornPhSpPt& born, double muR2) const override;

        std::vector<BornChannel> bornChannels() const override;
        
    private:        
        double bornAmp2(int flavour, const FourVector& pQ, const FourVector& pQbar, const FourVector& pLMinus, 
                const FourVector& pLPlus) const;

        double realAmp2_crossed( int flavour, const FourVector& qSlot, const FourVector& qbarSlot, const FourVector& gluonSlot, 
                const FourVector& pLMinus, const FourVector& pLPlus, double alphaS, double colorFactor) const;

        double realAmp2qqbar(const RealPhSpPt real, double muR2, bool useCMWAlphaS) const;
        double realAmp2gluonLeg1(const RealPhSpPt real, double muR2, bool useCMWAlphaS) const;
        double realAmp2gluonLeg2(const RealPhSpPt real, double muR2, bool useCMWAlphaS) const;

        double realOverBornQQbar(const RealPhSpPt& real, double muR2, bool useCMWALphaS) const;
        double realOverBornGQbar(const RealPhSpPt& real, double muR2, bool useCMWALphaS) const;
        double realOverBornQG(const RealPhSpPt& real, double muR2, bool useCMWALphaS) const;
    };

} // namespace powheg_dy
