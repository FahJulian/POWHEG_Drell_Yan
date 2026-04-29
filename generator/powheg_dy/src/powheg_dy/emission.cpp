#include "emission.h"

#include "powheg_dy/rand.h"
#include "powheg_dy/process.h"
#include "powheg_dy/born_event.h"

#include <math.h>

namespace powheg_dy
{
    namespace 
    {
        static constexpr double __DELTA = 1.0e-6;
        static constexpr double __T_CUTOFF = 2;   // GeV^2
        static constexpr double __B = 1.5;

        // The integral over the z distribution
        double __computeZIntegral(double zMin, double zMax)
        {
            return std::log((1 - zMin) / (1 - zMax));
        }

        // Sample a random Sudakov factor and compute t from it
        double __sampleSudakovAndComputeT(double tMax, double normZ)
        {
            double sudakov = rand(0.0, 1.0);
            return tMax * exp(log(sudakov) / __B / normZ);
        }

        // Sample a random z from the appropriate probability distribution
        double __sampleZ(double zMin, double normZ)
        {
            double u = rand(0.0, 1.0);
            return 1.0 - (1.0 - zMin) * exp(-u * normZ);
        }

        // Compute the franction of events with the given kinematics to accept
        double __computeAcceptanceRatio(double t, double z, double x, int partonId, 
            const std::unique_ptr<LHAPDF::PDF>& pdf)
        {
            double pdfRatio = pdf->xfxQ2(partonId, x / z, t) / pdf->xfxQ2(partonId, x, t);
            double prefactor = Physics::alphaSOneLoop(t, 5) / 2 / Math::PI * Physics::C_F / __B;
            
            return prefactor * (1 + z * z) * pdfRatio;
        }
    }
    
    Emission Emission::generateFirstEmission(const Process& process, const BornEvent& bornEvent)
    {
        // For each leg, sample one emission
        Emission emission1 = _generateEmissionOnLeg(process, bornEvent, 1);
        Emission emission2 = _generateEmissionOnLeg(process, bornEvent, 2);

        // Choose the emission on the leg with the largest pT^2
        double weight1 = emission1.isRejected() ? 0.0 : emission1.m_t;
        double weight2 = emission2.isRejected() ? 0.0 : emission2.m_t;

        return weight1 > weight2 ? emission1 : emission2;
    }

    Emission Emission::_generateEmissionOnLeg(const Process& process, const BornEvent& bornEvent, int leg)
    {
        Emission emission(process, bornEvent, leg);

        // For now don't generate an event for charm or bottom bosons since for those B is huge
        if (std::abs(bornEvent.getPartonId()) > 3) 
            return emission.reject();

        double x = leg == 1 ? bornEvent.getX1() : bornEvent.getX2();

        // Determine the kinematic bounds on z: The upper bound is t-dependent, so we use no 
        // upper bound at first and then veto if the bound is exceeded
        double zMin = x;
        double zMax = 1.0 - __DELTA;

        // Compute the normalization of the z-distribution, i.e. the integral over it
        double normZ = __computeZIntegral(zMin, zMax);

        double tMax = bornEvent.getS();
        while (tMax > __T_CUTOFF)
        {
            // Sample a Sudakov factor and compute t from it, check if it is above the cutoff
            emission.m_t = __sampleSudakovAndComputeT(tMax, normZ);
            if (emission.m_t < __T_CUTOFF)
                return emission.reject();

            // Sample z from the appropriate probability distribution
            emission.m_z = __sampleZ(zMin, normZ);

            // Sample phi uniformly
            emission.m_phi = rand(0, 2.0 * Math::PI);

            // Compute the franction of events with the given kinematics to accept
            double rAcc = __computeAcceptanceRatio(emission.m_t, emission.m_z, x, 
                leg == 1 ? bornEvent.getPartonId() : -bornEvent.getPartonId(), process.getPdfs());
            ASSERT (rAcc <= 1.0, "Acceptance ratio is greater than one: Upper bound B too small");

            // If the kinematics are allowed, accept the emission with the ratio just computed, 
            // else try another emission at the lower scale
            double r = emission.m_t / bornEvent.getS();
            double zMaxPhys = std::pow(std::sqrt(1 + r) - std::sqrt(r), 2);
            if (emission.m_z < zMaxPhys && rand(0.0, 1.0) < rAcc)
                return emission;
            else 
                tMax = emission.m_t;
        }

        // No emission above cutoff could be generated
        return emission.reject();
    }

} // namespace powheg_dy
