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
        virtual inline double wMass() const { return 80.398; }
        virtual inline double zWidth() const { return 2.4952; }
        virtual inline double ALPHA() const { return 1.0/128.89; }
        virtual inline double NC() const { return 3.0; }
        virtual inline double GEV2_TO_PB() const { return 0.389379338e9; }
        virtual inline double S_W_SQ() const { return 1.0 - wMass()*wMass()/zMass()/zMass(); }
        virtual inline double C_W_SQ() const { return 1.0 - S_W_SQ(); }
        virtual inline double M_Z() const { return 91.1876; }   
        virtual inline double GAMMA_Z() const { return 2.4952; }    
        virtual inline double KAPPA() const { return 1.0 / (4.0 * S_W_SQ() * C_W_SQ()); }
        virtual inline double LAMBDA_SQ_QCD() const { return 0.22626125081841550*0.22626125081841550; }
        virtual inline double C_F() const { return 4.0 / 3.0; }
        virtual inline double C_A() const { return 3.0; }
        virtual inline double T_F() const { return 0.5; }
        virtual inline double CMW_CHARM_THRESHOLDSQ() const { return 1.5*1.5; }
        virtual inline double CMW_BOTTOM_THRESHOLDSQ() const { return 5.0*5.0; }

        // charged lepton axial and vector couplings
        virtual inline double A_L() const { return -0.5; }
        virtual inline double V_L() const { return -0.5 + 2.0 * S_W_SQ(); }

        virtual inline double alphaSFromPdf(double qSq) const
        {
            return m_pdfs->alphasQ2(qSq);
        }

        inline double beta0(int nF) const
        {
            return (11.0 * C_A() - 4.0 * T_F() * nF) / (12.0 * PI);
        }

        inline double alphaS0(double qSq, int nF, double lambdaSq=-1.0) const
        {
            if (lambdaSq == -1.0)
                lambdaSq = LAMBDA_SQ_QCD();

            return 1.0 / (beta0(nF) * std::log(qSq / lambdaSq));
        }

        inline int activeFlavours(double qSq) const
        {
            if (qSq < CMW_CHARM_THRESHOLDSQ()) 
                return 3;
            else if (qSq < CMW_BOTTOM_THRESHOLDSQ()) 
                return 4;
            else 
                return 5;
        }

        inline double beta0Powheg(int nf) const
        {
            return (33.0 - 2.0 * nf) / (12.0 * PI);
        }

        inline double betaPrimePowheg(int nf) const
        {
            return (153.0 - 19.0 * nf) / (2.0 * PI * (33.0 - 2.0 * nf));
        }

        inline double alphaRawPowheg(double qSq, double lambda5Sq, int nf) const
        {
            const double L = std::log(qSq / lambda5Sq);
            assert(L > 0.0);

            const double b  = beta0Powheg(nf);
            const double bp = betaPrimePowheg(nf);

            return 1.0 / (b * L) - bp * std::log(L) / ((b * L) * (b * L));
        }

        inline double alphaSPowheg(double qSq) const
        {
            const double lambda5Sq = LAMBDA_SQ_QCD();

            const double mc2 = CMW_CHARM_THRESHOLDSQ();
            const double mb2 = CMW_BOTTOM_THRESHOLDSQ();

            const double a5_mb_raw = alphaRawPowheg(mb2, lambda5Sq, 5);
            const double a4_mb_raw = alphaRawPowheg(mb2, lambda5Sq, 4);

            const double c45 = 1.0 / a5_mb_raw - 1.0 / a4_mb_raw;

            const double a4_mc_raw = alphaRawPowheg(mc2, lambda5Sq, 4);
            const double a3_mc_raw = alphaRawPowheg(mc2, lambda5Sq, 3);

            const double c35 = 1.0 / a4_mc_raw - 1.0 / a3_mc_raw + c45;

            const int nf = activeFlavours(qSq);

            if (nf == 5)
            {
                return alphaRawPowheg(qSq, lambda5Sq, 5);
            }
            else if (nf == 4)
            {
                const double a4_raw = alphaRawPowheg(qSq, lambda5Sq, 4);
                return 1.0 / (1.0 / a4_raw + c45);
            }
            else
            {
                const double a3_raw = alphaRawPowheg(qSq, lambda5Sq, 3);
                return 1.0 / (1.0 / a3_raw + c35);
            }
        }

        inline double alphaSCMW(double qSq) const
        {
            const int nF = activeFlavours(qSq);
            // const double alphaS = alphaSFromPdf(qSq);
            const double alphaS = alphaSPowheg(qSq);

            const double K = (67.0 / 18.0 - PI * PI / 6.0) * C_A() - 5.0 / 9.0 * nF;

            return alphaS * (1.0 + alphaS / (2.0 * PI) * K);
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
