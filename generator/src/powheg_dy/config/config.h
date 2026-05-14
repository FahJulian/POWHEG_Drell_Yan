#pragma once

#include "powheg_dy/base.h"

namespace powheg_dy
{
    struct Config
    {
        ////////// mandatory params ///////////
        double SQRT_S;

        double M_MIN;
        double M_MAX;

        double PT_SQ_CUTOFF;

        size_t N_ACCEPTED_EVENTS;
        size_t N_TRIAL_EVENTS;

        std::string PDF_NAME;

        ////////// other independent params ///////////
        bool NO_EMISSIONS = false;
        bool BTILDE_BORNONLY = false;

        bool ALPHA_S_FROM_PDF = false;

        double BORN_VETO_WEIGHT = -1.0;
        
        /////////// dependent params ////////////
        double S;

        ///////// independent constants /////////
        double M_Z = 91.1876;
        double M_W = 80.398;
        double GAMMA_Z = 2.4952;
        double GAMMA_W = 2.141;

        double ALPHA_EW = 1.0 / 128.89;

        double N_C = 3.0;
        double C_F = 4.0 / 3.0;
        double C_A = 3.0;
        double T_F = 0.5;

        double CMW_CHARM_TRSHLD_SQ = 1.5 * 1.5;
        double CMW_BOTTOM_TRSHLD_SQ = 5.0 * 5.0;

        ////////// dependent constants //////////
        double S_W_SQ;
        double C_W_SQ;
        double S_W;
        double C_W;
        
        double E_SQ;
        
        ////////////// set on init ///////////////
        double LAMBDA_MSB_5_SQ;
        std::unique_ptr<LHAPDF::PDF> PDF;

        double beta0(int nF) const
        {
            return (11.0 * C_A - 4.0 * T_F * nF) / (12.0 * PI);
        }
    
        double alphaS0(double qSq, int nF) const;
        double alphaS0customLambda(double qSq, int nF, double lambdaSq) const;
        double alphaS(double qSq) const;
        double alphaSCMW(double qSq) const;

    private:
        void initLHAPDF();
        void setDependentParams();
        void extractLambdaFromPdf();

        friend class Process;
    };


} // namespace powheg_dy
