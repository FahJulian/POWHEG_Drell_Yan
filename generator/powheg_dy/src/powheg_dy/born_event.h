#pragma once

#include "powheg_dy/base.h"

#include <tuple>
#include <vector>

namespace powheg_dy
{
    class Process;

    class BornEvent
    {
    public:
        BornEvent() = delete;
        BornEvent(const Process& process)
            : m_process(process)
        {
        }

        void sampleKinematics();
        void computeWeightAndSampleParton();

        double getMBoson() const { return m_mBoson; }
        double getYBoson() const { return m_yBoson; }
        double getS() const { return m_mBoson * m_mBoson; }
        double getX1() const { return m_x1; }
        double getX2() const { return m_x2; }
        double getPhi() const { return m_phi; }
        double getCosTh() const { return m_cosTh; }
        double getDSigma() const { return m_dSigma; }
        int getPartonId() const { return m_partonId; }
    
    private:
        double _computeInverseSamplingDensity();
        std::vector<std::tuple<int, double>> _computePartonChannelContributions();

    private:
        const Process& m_process;
        int m_partonId = 0;
        double m_mBoson = 0.0;
        double m_yBoson = 0.0;
        double m_phi = 0.0;
        double m_cosTh = 0.0;   // The angle between the lepton and the quark, not between the lepton and leg 1!!
        double m_x1, m_x2;
        double m_dSigma = 0.0;
    };

} // namespace powheg_dy
