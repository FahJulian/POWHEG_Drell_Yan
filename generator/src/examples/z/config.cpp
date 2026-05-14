#include "drell_yan.h"

namespace powheg_dy
{
    void DrellYanProcess::initConfig(ConfigParser& parser) const
    {
        powheg_assert(parser.extract("M_MIN", m_config.M_MIN), "No value for M_MIN found.");
        powheg_assert(parser.extract("M_MAX", m_config.M_MAX), "No value for M_MAX found.");

        parser.extract("M_Z", m_config.M_Z);
        parser.extract("M_W", m_config.M_W);
        parser.extract("GAMMA_Z", m_config.GAMMA_Z);
        parser.extract("GAMMA_W", m_config.GAMMA_W);

        parser.extract("ALPHA_EW", m_config.ALPHA_EW);

        m_config.M_Z_SQ = m_config.M_Z * m_config.M_Z;
        m_config.M_W_SQ = m_config.M_W * m_config.M_W;

        m_config.S_W_SQ = 1.0 - m_config.M_W_SQ / m_config.M_Z_SQ;
        m_config.C_W_SQ = 1.0 - m_config.S_W_SQ;

        m_config.S_W = std::sqrt(m_config.S_W_SQ);
        m_config.C_W = std::sqrt(m_config.C_W_SQ);

        m_config.E_SQ = 4.0 * PI * m_config.ALPHA_EW;
    }

    ZHelicityCouplings DrellYanConfig::zCouplings(const bool upType, const double charge) const
    {
        int sign = upType ? 1 : -1;

        const double v = 0.5 * sign - 2.0 * charge * S_W_SQ;
        const double a = 0.5 * sign;

        ZHelicityCouplings out;
        out.minus = v + a;      // Fortran Zcoup(-1)
        out.plus  = v - a;      // Fortran Zcoup(+1)
        return out;
    }

    std::complex<double> DrellYanConfig::photonPropagator(const double q2) const
    {
        return 1.0 / q2;
    }

    std::complex<double> DrellYanConfig::zPropagator(const double q2) const
    {
        // no running width
        return 1.0 / std::complex<double>(q2 - M_Z_SQ, M_Z * GAMMA_Z);
    }
    
} // namespace powheg_dy
