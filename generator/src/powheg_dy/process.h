#pragma once

#include "powheg_dy/base.h"
#include "powheg_dy/event.h"
#include "powheg_dy/config.h"
#include "powheg_dy/matrix_elements.h"
#include "powheg_dy/bbar/bbar_integrator.h"
#include "powheg_dy/emission/emission_generator.h"
#include "powheg_dy/phase_space/born_phase_space.h"
#include "powheg_dy/phase_space/real_phase_space.h"
#include "powheg_dy/born_event/born_event_generator.h"

namespace powheg_dy
{
    class Process
    {
    public:
        ~Process();

        virtual double bornAmp2(const BornPhSpPt& born) const = 0;
        virtual double virtAmp2(const BornPhSpPt& born, const double amp2Born, const double muR2) const = 0;
        virtual double realAmp2(const RealPhSpPt& real, const RealChannel& channel, const double alphaS) const = 0;

        virtual std::vector<BornChannel> bornChannels() const = 0;
        virtual std::vector<RealChannel> realChannels(const BornChannel& bornChannel) const = 0;

        // TODO: Move born sampling here
        
        // OLD STUFF
        virtual RealOverBornContributions realOverBornContributions(
            const RealPhSpPt& real, double muF2, double muR2, bool useCMWALphaS) const = 0;

        virtual double bornContribution(const BornPhSpPt& born) const = 0;
        virtual double virtualContribution(const BornPhSpPt& born, double muR2) const = 0;

        void init(const std::string& configPath);
        void run();
        void writeToFile(const std::string& filePath) const;

        double getSigma() const { return m_bbarIntegrator->getTotalCrossSection(); }
        auto getEvents() const { return m_events; }

    private:
        void clear();
        void analyse();
        void generateEvents();

    private:
        bool m_initialized = false;

        std::shared_ptr<BornPhaseSpace> m_bornPhSp;
        std::shared_ptr<FKSRealPhaseSpace> m_realPhSp;
        std::shared_ptr<BBarIntegrator> m_bbarIntegrator;
        std::shared_ptr<BornEventGenerator> m_bornGenerator;
        std::shared_ptr<EmissionGenerator> m_emissionGenerator;

        std::vector<Event> m_events;

    protected:
        Config m_config;
    };  

} // namespace powheg_dy
