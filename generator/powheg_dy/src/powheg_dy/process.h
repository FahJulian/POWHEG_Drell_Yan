#pragma once

#include "powheg_dy/math/math.h"
#include "powheg_dy/event.h"
#include "powheg_dy/born_event/born_event_generator.h"
#include "powheg_dy/phase_space/born_phase_space.h"
#include "powheg_dy/phase_space/real_phase_space.h"
#include "powheg_dy/emission/emission_generator.h"

#include <memory>
#include <vector>
#include <string>

#include <LHAPDF/LHAPDF.h>

namespace powheg_dy
{
    class Process
    {
    public:
        virtual ~Process() = default;

        void init(const std::string& pdfDataLocation, const std::string& pdfSet);
        void run();
        void writeToFile(const std::string& filePath) const;

        virtual inline int nBornLegs() const = 0;
        virtual inline int nRealLegs() const = 0;
        virtual inline bool bornOnly() const = 0;
        virtual inline double mMin() const = 0;
        virtual inline double mMax() const = 0;
        virtual inline double sqrtS() const = 0;
        virtual inline double S() const = 0;
        virtual inline double pt2Cutoff() const = 0;

        virtual inline double zMass() const { return 91.1876; }
        virtual inline double zWidth() const { return 2.4952; }
        virtual inline double ALPHA() const { return 1.0 / 137.035999084; }
        virtual inline double NC() const { return 3.0; }
        virtual inline double GEV2_TO_PB() const { return 0.389379338e9; }
        virtual inline double S_W_SQ() const { return 0.23126; }
        virtual inline double C_W_SQ() const { return 1.0 - S_W_SQ(); }
        virtual inline double M_Z() const { return 91.1876; }   
        virtual inline double GAMMA_Z() const { return 2.4952; }    
        virtual inline double KAPPA() const { return 1.0 / (4.0 * S_W_SQ() * C_W_SQ()); }
        virtual inline double LAMBDA_SQ_QCD() const { return 0.2*0.2; }
        virtual inline double C_F() const { return 4.0 / 3.0; }
        virtual inline double C_A() const { return 3.0; }
        virtual inline double T_F() const { return 0.5; }

        // charged lepton axial and vector couplings
        virtual inline double A_L() const { return -0.5; }
        virtual inline double V_L() const { return -0.5 + 2.0 * S_W_SQ(); }

        virtual inline double alphaSOneLoop(double qSq, int nF) const
        {
            assert(qSq > LAMBDA_SQ_QCD());

            double beta0 = 11.0 - 2.0 / 3.0 * nF;
            return 4.0 * PI / beta0 / std::log(qSq / LAMBDA_SQ_QCD());
        }
        
        double getSigma() const { return m_totalCrossSection; }
        auto getEvents() const { return m_events; }
        const std::unique_ptr<LHAPDF::PDF>& getPdfs() const { return m_pdfs; }

    private:
        void _clear();
        void _determineMaxWeight();
        void _computeTotalCrossSection();

    private:
        std::unique_ptr<BornPhaseSpace> m_bornPhSp;
        std::unique_ptr<FKSRealPhaseSpace> m_realPhSp;
        std::unique_ptr<BornEventGenerator> m_bornGenerator;
        std::unique_ptr<EmissionGenerator> m_emissionGenerator;
        int m_nEventTrials = 0;
        double m_maxWeight = 0.0;
        double m_totalCrossSection = 0.0;   // pb

        std::unique_ptr<LHAPDF::PDF> m_pdfs;
        std::vector<Event> m_events;
    };  

} // namespace powheg_dy
