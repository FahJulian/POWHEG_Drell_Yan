#include "emission.h"

#include "powheg_dy/math.h"
#include "powheg_dy/rand.h"
#include "powheg_dy/process.h"
#include "powheg_dy/born_event.h"


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
            return log((1 - zMin) / (1 - zMax));
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
        double __computeAcceptanceRatio(const Process& process, double t, double z, double x, int partonId, 
            const std::unique_ptr<LHAPDF::PDF>& pdf)
        {
            double pdfRatio = pdf->xfxQ2(partonId, x / z, t) / pdf->xfxQ2(partonId, x, t);
            double prefactor = process.alphaSOneLoop(t, 5) / 2 / PI * process.C_F() / __B;
            
            return prefactor * (1 + z * z) * pdfRatio;
        }
    }
    
    Emission EmissionGenerator::generateEmission(const BornPhaseSpacePt& point) const
    {
        // For each leg, sample one emission
        Emission emission1 = _generateEmissionOnLeg(point, 1);
        Emission emission2 = _generateEmissionOnLeg(point, 2);

        // Choose the emission on the leg with the largest pT^2
        double weight1 = emission1.rejected ? 0.0 : emission1.t;
        double weight2 = emission2.rejected ? 0.0 : emission2.t;

        return weight1 > weight2 ? emission1 : emission2;
    }

    Emission EmissionGenerator::_generateEmissionOnLeg(const BornPhaseSpacePt& point, int leg) const
    {
        Emission emission;
        emission.leg = leg;

        // For now don't generate an event for charm or bottom bosons since for those B is huge
        if (point.channel.flavour > 3) 
            return emission.reject();

        double x = leg == 1 ? point.x1Bar : point.x2Bar;

        // Determine the kinematic bounds on z: The upper bound is t-dependent, so we use no 
        // upper bound at first and then veto if the bound is exceeded
        double zMin = x;
        double zMax = 1.0 - __DELTA;

        // Compute the normalization of the z-distribution, i.e. the integral over it
        double normZ = __computeZIntegral(zMin, zMax);

        double tMax = point.sHat;
        while (tMax > __T_CUTOFF)
        {
            // Sample a Sudakov factor and compute t from it, check if it is above the cutoff
            emission.t = __sampleSudakovAndComputeT(tMax, normZ);
            if (emission.t < __T_CUTOFF)
                return emission.reject();

            // Sample z from the appropriate probability distribution
            emission.z = __sampleZ(zMin, normZ);

            // Sample phi uniformly
            emission.phi = rand(0, 2.0 * PI);

            // Compute the franction of events with the given kinematics to accept
            double rAcc = __computeAcceptanceRatio(m_process, emission.t, emission.z, x, 
                leg == 1 ? point.channel.id1 : point.channel.id2, m_process.getPdfs());
            assert(rAcc <= 1.0);    // Acceptance ratio is greater than one: Upper bound B too small

            // If the kinematics are allowed, accept the emission with the ratio just computed, 
            // else try another emission at the lower scale
            double r = emission.t / point.sHat;
            double zMaxPhys = pow(sqrt(1 + r) - sqrt(r), 2);
            if (emission.z < zMaxPhys && rand(0.0, 1.0) < rAcc)
                return emission;
            else 
                tMax = emission.t;
        }

        // No emission above cutoff could be generated
        return emission.reject();
    }

} // namespace powheg_dy
