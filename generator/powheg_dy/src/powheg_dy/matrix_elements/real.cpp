#include "matrix_elements.h"

#include "powheg_dy/math.h"
#include "powheg_dy/process.h"

namespace powheg_dy
{
    namespace 
    {
        constexpr int __NF = 5;
        constexpr double __EPS = 1.0e-14;

        double __Fqqbar(double xi, double y)
        {
            const double oneMinusY2 = max(1.0 - y * y, __EPS);
            const double xi2 = max(xi * xi, __EPS);

            // The bracket in eq. (7.200)
            return 2.0 * (1.0 + y * y) / oneMinusY2
                + 8.0 * (1.0 - xi) / (xi2 * oneMinusY2);
        }

        double __Fgqbar(double xi, double y)
        {
            const double denom = max(xi * (1.0 - y), __EPS);

            // The negative of the bracket in eq. (7.200)
            return 2.0 * (1.0 - xi * (1.0 - xi) * (1.0 + y)) / denom 
                + 0.5 * xi * (1.0 - y);
        }

    } // namespace

    double MatrixElements::realOverBornQQbar(const Process& process, const RealPhSpPt& real, double muR2)
    {
        const double gs2 = 4.0 * PI * process.alphaSOneLoop(muR2, __NF);

        // eq. (7.200), divided by the born factor g^2/4N_C
        return 4.0 * process.C_F() * gs2 / real.sHatReal * __Fqqbar(real.rad.xi, real.rad.y);
    }

    double MatrixElements::realOverBornGQbar(const Process& process, const RealPhSpPt& real, double muR2)
    {
        const double gs2 = 4.0 * PI * process.alphaSOneLoop(muR2, __NF);

        // eq. (7.201), divided by the born factor g^2/4N_C
        return 4.0 * process.T_F() * gs2 / real.sHatReal * __Fgqbar(real.rad.xi, real.rad.y);
    }

    double MatrixElements::realOverBornQG(const Process& process, const RealPhSpPt& real, double muR2)
    {
        const double gs2 = 4.0 * PI * process.alphaSOneLoop(muR2, __NF);

        // The same as realOverBornGQbar, but with y \to -y
        return 4.0 * process.T_F() * gs2 / real.sHatReal * __Fgqbar(real.rad.xi, -real.rad.y);
    }

    MatrixElements::RealOverBornContributions MatrixElements::realOverBornContributions(
        const Process& process,
        const RealPhSpPt& real,
        double muF2,
        double muR2
    )
    {
        RealOverBornContributions out;

        const BornPhSpPt& born = real.underlyingBorn;
        const int id1 = born.channel.id1;
        const int id2 = born.channel.id2;

        const double fBorn1 = process.getPdfs()->xfxQ2(id1, born.x1Bar, muF2) / born.x1Bar;
        const double fBorn2 = process.getPdfs()->xfxQ2(id2, born.x2Bar, muF2) / born.x2Bar;
        const double bornLuminosity = fBorn1 * fBorn2;

        const double fReal1Q  = process.getPdfs()->xfxQ2(id1, real.x1, muF2) / real.x1; 
        const double fReal2QB = process.getPdfs()->xfxQ2(id2, real.x2, muF2) / real.x2; 

        const double fReal1G  = process.getPdfs()->xfxQ2(21, real.x1, muF2) / real.x1; 
        const double fReal2G  = process.getPdfs()->xfxQ2(21, real.x2, muF2) / real.x2; 

        const double lumQQbar = fReal1Q * fReal2QB;
        const double lumGQbar = fReal1G * fReal2QB;
        const double lumQG    = fReal1Q * fReal2G;

        out.qqbar = lumQQbar / bornLuminosity * realOverBornQQbar(process, real, muR2);
        out.gqbar = lumGQbar / bornLuminosity * realOverBornGQbar(process, real, muR2);
        out.qg = lumQG / bornLuminosity * realOverBornQG(process, real, muR2);

        return out;
    }

    double MatrixElements::realOverBorn(
        const Process& process,
        const RealPhSpPt& real,
        double muF2,
        double muR2
    )
    {
        return realOverBornContributions(process, real, muF2, muR2).total();
    }

    double MatrixElements::realSudakovDensity(
        const Process& process,
        const RealPhSpPt& real,
        double muF2,
        double muR2
    )
    {
        return real.radJacobian * realOverBorn(process, real, muF2, muR2);
    }

} // namespace powheg_dy