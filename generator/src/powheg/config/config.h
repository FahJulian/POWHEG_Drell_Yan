#pragma once

#include "powheg/base.h"

namespace powheg
{
    struct Config
    {
        virtual ~Config() = default;

        ////////// mandatory process params ///////////
        double SQRT_S;

        double PT_SQ_CUTOFF;

        size_t N_ACCEPTED_EVENTS;
        size_t N_TRIAL_EVENTS;

        std::string PDF_NAME;

        ////////// other independent process params ///////////
        bool NO_EMISSIONS = false;
        bool BTILDE_BORNONLY = false;

        bool BTILDE_USE_BORN = true;
        bool BTILDE_USE_VIRTUAL = true;
        bool BTILDE_USE_REAL = true;
        bool BTILDE_USE_COUNTERTERMS = true;
        bool BTILDE_USE_COLL_REMNANTS = true;

        bool ALPHA_S_FROM_PDF = false;

        double BORN_VETO_WEIGHT = -1.0;
        
        ////////////// dependent process params ///////////////
        double S;

        //////////// independent physics constants ////////////
        double N_C = 3.0;
        double C_F = 4.0 / 3.0;
        double C_A = 3.0;
        double T_F = 0.5;

        double CMW_CHARM_TRSHLD = 1.5;
        double CMW_BOTTOM_TRSHLD = 5.0;

        ///////////// dependent physics constants /////////////
        double CMW_CHARM_TRSHLD_SQ;
        double CMW_BOTTOM_TRSHLD_SQ;
        
        ///////////////////// set on init /////////////////////
        double LAMBDA_MSB_5_SQ;
        std::unique_ptr<LHAPDF::PDF> PDF;

        double beta0(int nF) const;
        double alphaS0(double qSq, int nF) const;
        double alphaS0customLambda(double qSq, int nF, double lambdaSq) const;
        double alphaS(double qSq) const;
        double alphaSCMW(double qSq) const;

    private:
        void initPdf();
        void setDependentParams();
        void extractLambdaFromPdf();

        friend class BaseProcess;
    };


} // namespace powheg
