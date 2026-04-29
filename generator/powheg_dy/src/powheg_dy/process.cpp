#include "process.h"

#include "powheg_dy/rand.h"
#include "powheg_dy/file.h"
#include "powheg_dy/emission.h"
#include "powheg_dy/born_event.h"
#include "powheg_dy/les_houches_serializer.h"

#include <iostream>

namespace powheg_dy
{
    namespace 
    {
        static constexpr int __N_TRIAL_EVENTS = 1e5;
        static constexpr int __N_ACCEPTED_EVENTS = 1e4;
        static constexpr double __SECURITY_FACTOR = 1.1;

    } // namespace
    
    void Process::init(const std::string& pdfDataLocation, const std::string& pdfSet)
    {
        LHAPDF::setPaths(pdfDataLocation);
        m_pdfs = std::unique_ptr<LHAPDF::PDF>(LHAPDF::mkPDF(pdfSet, 0));
    }

    void Process::run()
    {
        _clear();

        _determineMaxWeight();
        
        while (m_events.size() < __N_ACCEPTED_EVENTS)
        {   
            m_nEventTrials++;

            BornEvent bornEvent(*this);
            bornEvent.sampleKinematics();
            bornEvent.computeWeightAndSampleParton();

            ASSERT(bornEvent.getDSigma() <= m_maxDSigma, "Born weight exceeds max weight");
            
            // Unweight: Accept the event at the rate of it's weight over the reference weight
            double u = rand(0.0, 1.0);
            if (u < bornEvent.getDSigma() / m_maxDSigma)
            {
                auto emission = Emission::generateFirstEmission(*this, bornEvent);

                Event event(*this, bornEvent, emission);
                event.reconstructMomenta();
                
                m_events.push_back(event);
            }
        }

        std::cout << "Acceptance ratio: " << double(__N_ACCEPTED_EVENTS) / m_nEventTrials << std::endl;
        
        _computeTotalCrossSection();
        
        std::cout << "Total cross section: " << m_totalCrossSection << " pb." << std::endl;
    }

    void Process::writeToFile(const std::string& filePath) const
    {
        std::string fileContent;
            
        for (const auto& event : m_events)
            fileContent.append(event.toString() + '\n');
        
        File file = File(filePath);
        file.write(fileContent);

        LesHouchesSerializer(*this).serialize(filePath);
    }

    void Process::_clear()
    {
        m_events.clear();
        m_events.reserve(__N_ACCEPTED_EVENTS);
        m_nEventTrials = 0;
        m_totalCrossSection = 0.0;
        m_maxDSigma = 0.0;
    }

    void Process::_determineMaxWeight()
    {
        double max_dSigma = 0.0;
        
        for (int i = 0; i < __N_TRIAL_EVENTS; i++)
        {   
            BornEvent bornEvent(*this);
            bornEvent.sampleKinematics();
            bornEvent.computeWeightAndSampleParton();

            if (bornEvent.getDSigma() > max_dSigma)
                max_dSigma = bornEvent.getDSigma();
        }

        m_maxDSigma = __SECURITY_FACTOR * max_dSigma;
    }

    void Process::_computeTotalCrossSection()
    {
        m_totalCrossSection = m_maxDSigma * __N_ACCEPTED_EVENTS / m_nEventTrials * Physics::GEV2_TO_PB;
    }

} // namespace powheg_dy
