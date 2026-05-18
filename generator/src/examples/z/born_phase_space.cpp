#include "drell_yan.h"

#include "powheg_dy/math/rand.h"

namespace powheg
{
    BornPhSpPt DrellYanProcess::sampleBorn() const
    {
        const auto [mSq, jacobianMSq] = sampleMSq();
        const auto [yBoson, jacobianYBoson] = sampleYBoson(mSq);
        const auto [cosTh, jacobianCosTh] = sampleCosTh();
        const auto [phi, jacobianPhi] = samplePhi();

        const double mBoson = std::sqrt(mSq);
        const double sinTh = std::sqrt(1.0 - cosTh * cosTh);
        const double sinPhi = std::sin(phi);
        const double cosPhi = std::cos(phi);

        BornPhSpPt born;

        born.sHat = mSq;
        born.x1Bar = mBoson / m_config.SQRT_S * std::exp(yBoson);
        born.x2Bar = mBoson / m_config.SQRT_S * std::exp(-yBoson);

        born.jacobian = 1.0 / (32.0 * PI * PI)      // d\Phi = 1/32\pi^2 d\cos\theta d\phi
            * 1.0 / m_config.S                      // change of variables from (x1, x2) to (sHat, yB)
            * jacobianMSq * jacobianYBoson          // sampling jacobians
            * jacobianCosTh * jacobianPhi
        ;   

        born.p1Bar = {
            0.5 * born.x1Bar * m_config.SQRT_S,
            0.0,
            0.0,
            0.5 * born.x1Bar * m_config.SQRT_S
        };

        born.p2Bar = {
            0.5 * born.x2Bar * m_config.SQRT_S,
            0.0,
            0.0,
            -0.5 * born.x2Bar * m_config.SQRT_S
        };

        const FourVector pLMinusRest = {
            mBoson / 2.0,
            mBoson / 2.0 * sinTh * sinPhi,
            mBoson / 2.0 * sinTh * cosPhi,
            mBoson / 2.0 * cosTh
        };

        const FourVector pLPlusRest = {
            mBoson / 2.0,
            -mBoson / 2.0 * sinTh * sinPhi,
            -mBoson / 2.0 * sinTh * cosPhi,
            -mBoson / 2.0 * cosTh
        };
        
        const FourVector pBoson = { mBoson * std::cosh(yBoson), 0.0, 0.0, mBoson * std::sinh(yBoson) };
        const ThreeVector betaBoson = pBoson.getBeta();

        born.pOut.push_back(pLMinusRest.boost(betaBoson));
        born.pOut.push_back(pLPlusRest. boost(betaBoson));
        
        return born;
    }

    std::pair<double, double> DrellYanProcess::sampleMSq() const
    {
        const double u = rand();

        const double m2Low = m_config.M_MIN * m_config.M_MIN;
        const double m2High = m_config.M_MAX * m_config.M_MAX;

        // Breit-Wigner-like sampling of m^2 around the Z peak
        const double zlow = std::atan((m2Low - m_config.M_Z_SQ) / (m_config.M_Z * m_config.GAMMA_Z));
        const double zhigh = std::atan((m2High - m_config.M_Z_SQ) / (m_config.M_Z * m_config.GAMMA_Z));

        const double z = zlow + (zhigh - zlow) * u;
        const double mSq = m_config.M_Z * m_config.GAMMA_Z * std::tan(z) + m_config.M_Z_SQ;

        // Jacobian dm^2 / du
        const double jacobianMSq = (zhigh - zlow) * m_config.M_Z * m_config.GAMMA_Z 
            / (std::cos(z) * std::cos(z));

        return { mSq, jacobianMSq };
    }

    std::pair<double, double> DrellYanProcess::sampleYBoson(const double sHat) const
    {
        const double u = rand();

        // Sample boson rapidity uniformly from the kinematically allowed range
        const double yBosonMax = 0.5 * std::log(m_config.S / sHat);
        const double yBoson = (2.0 * u - 1.0) * yBosonMax;

        const double jacobianYBoson = 2 * yBosonMax;

        return { yBoson, jacobianYBoson };
    }

    std::pair<double, double> DrellYanProcess::sampleCosTh() const
    {
        const double u = rand();

        // Sample cos(theta) from p(c) = 3(1 + c^2) / 8
        const double cosTh = 2.0 * std::sinh(std::asinh(4.0 * u - 2.0) / 3.0);
        const double jacobianCosTh = 8.0 / 3.0 / (1.0 + cosTh * cosTh);

        return { cosTh, jacobianCosTh };
    }

    std::pair<double, double> DrellYanProcess::samplePhi() const
    {
        const double u = rand();

        const double phi = 2.0 * PI * u;
        const double jacobianPhi = 2.0 * PI;

        return { phi, jacobianPhi };
    }

} // namespace powheg
