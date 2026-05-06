#pragma once

#include "powheg_dy/base.h"
#include "powheg_dy/event.h"
#include "powheg_dy/config.h"
#include "powheg_dy/emission/emission_generator.h"
#include "powheg_dy/phase_space/born_phase_space.h"
#include "powheg_dy/phase_space/real_phase_space.h"
#include "powheg_dy/born_event/born_event_generator.h"

#include <memory>
#include <vector>
#include <string>

#include <LHAPDF/LHAPDF.h>

namespace powheg_dy
{
    class Process
    {
    public:
        ~Process();
        
        void init(const std::string& pdfSet);
        void run();
        void writeToFile(const std::string& filePath) const;

        virtual inline int nBornLegs() const = 0;
        virtual inline int nRealLegs() const = 0;
        
        double getSigma() const { return m_totalCrossSection; }
        auto getEvents() const { return m_events; }

    private:
        void clear();
        void determineMaxWeight();
        void computeTotalCrossSection();

    private:
        std::shared_ptr<BornPhaseSpace> m_bornPhSp;
        std::shared_ptr<FKSRealPhaseSpace> m_realPhSp;
        std::shared_ptr<BornEventGenerator> m_bornGenerator;
        std::shared_ptr<EmissionGenerator> m_emissionGenerator;

        int m_nEventTrials = 0;
        double m_maxWeight = 0.0;
        double m_totalCrossSection = 0.0;   // pb

        std::vector<Event> m_events;

        Config m_config;
    };  

} // namespace powheg_dy
