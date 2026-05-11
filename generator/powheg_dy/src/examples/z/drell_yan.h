#pragma once

#include "powheg_dy/process.h"

namespace powheg_dy
{
    class DrellYanProcess : public Process
    {
    public:
        RealOverBornContributions realOverBornContributions(
            const RealPhSpPt& real, double muF2, double muR2, bool useCMWALphaS) const override;

        double born(const BornPhSpPt& born) const override;
        double virtualOverBorn(const BornPhSpPt& born, double muR2) const override;

        std::vector<BornChannel> bornChannels() const override;
        
    private:        
        double bornAmp2(int flavour, const FourVector& pQ, const FourVector& pQbar, const FourVector& pLMinus, 
                const FourVector& pLPlus) const;
        double realAmp2_crossed( int flavour, const FourVector& qSlot, const FourVector& qbarSlot, const FourVector& gluonSlot, 
                const FourVector& pLMinus, const FourVector& pLPlus, double alphaS, double colorFactor) const;

        double realOverBornQQbar(const RealPhSpPt& real, double muR2, bool useCMWALphaS) const;
        double realOverBornGQbar(const RealPhSpPt& real, double muR2, bool useCMWALphaS) const;
        double realOverBornQG(const RealPhSpPt& real, double muR2, bool useCMWALphaS) const;
    };

} // namespace powheg_dy
