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
        static constexpr int __N_TRIAL_EVENTS = 5e5;
        static constexpr int __N_ACCEPTED_EVENTS = 5e4;
        static constexpr double __SECURITY_FACTOR = 1.1;

    } // namespace
    
    void Process::init(const std::string& pdfDataLocation, const std::string& pdfSet)
    {
        LHAPDF::setPaths(pdfDataLocation);
        m_pdfs = std::unique_ptr<LHAPDF::PDF>(LHAPDF::mkPDF(pdfSet, 0));

        m_phaseSpaceSampler = std::make_unique<BornPhaseSpace>(*this);
        m_emissionGenerator = std::make_unique<EmissionGenerator>(*this);
        m_bornGenerator = std::make_unique<BornEventGenerator>(*this);
        m_eventHandler = std::make_unique<EventHandler>(*this);
    }

    void Process::run()
    {
        _clear();

        _determineMaxWeight();
        
        while (m_events.size() < __N_ACCEPTED_EVENTS)
        {   
            m_nEventTrials++;
            
            double rands[3] = { rand(), rand(), rand() };
            BornPhaseSpacePt point = m_phaseSpaceSampler->samplePoint(rands);
            m_bornGenerator->computeWeightAndSampleChannel(point);
            m_phaseSpaceSampler->reconstructMomenta(point);
            
            assert(point.weight <= m_maxWeight);
            
            // Unweight: Accept the event at the rate of it's weight over the reference weight
            double u = rand(0.0, 1.0);
            if (u < point.weight / m_maxWeight)
            {   
                if (bornOnly())
                {
                    Event event = m_eventHandler->reconstructEvent(point, Emission().reject());

                    m_events.push_back(event);
                }
                else
                {
                    Emission emission = m_emissionGenerator->generateEmission(point);
                    Event event = m_eventHandler->reconstructEvent(point, emission);

                    m_events.push_back(event);
                }
            }
        }

        std::cout << "Acceptance ratio: " << double(__N_ACCEPTED_EVENTS) / m_nEventTrials << std::endl;
        
        _computeTotalCrossSection();
        
        std::cout << "Total cross section: " << m_totalCrossSection << " pb." << std::endl;
    }

    std::string toString(const Event& event)
    {
        FourVector sum = event.p1Out + event.p2Out;
        double pT = sqrt(sum.pX*sum.pX + sum.pY*sum.pY);

        return std::to_string(event.point.mB) + ", "
            + std::to_string(event.point.cosTh) + ", "
            + std::to_string(pT);
    }

    void Process::writeToFile(const std::string& filePath) const
    {
        // std::string fileContent;
            
        // for (const auto& event : m_events)
        //     fileContent.append(toString(event) + '\n');
        
        // File file = File(filePath);
        // file.write(fileContent);

        LesHouchesSerializer(*this).serialize(filePath);
    }

    void Process::_clear()
    {
        m_events.clear();
        m_events.reserve(__N_ACCEPTED_EVENTS);
        m_nEventTrials = 0;
        m_totalCrossSection = 0.0;
        m_maxWeight = 0.0;
    }

    void Process::_determineMaxWeight()
    {
        double max_dSigma = 0.0;
        
        for (int i = 0; i < __N_TRIAL_EVENTS; i++)
        {   
            double rands[3] = { rand(), rand(), rand() };
            BornPhaseSpacePt point = m_phaseSpaceSampler->samplePoint(rands);
            m_bornGenerator->computeWeightAndSampleChannel(point);

            if (point.weight > max_dSigma)
                max_dSigma = point.weight;
        }

        m_maxWeight = __SECURITY_FACTOR * max_dSigma;
    }

    void Process::_computeTotalCrossSection()
    {
        m_totalCrossSection = m_maxWeight * __N_ACCEPTED_EVENTS / m_nEventTrials * GEV2_TO_PB();
    }

} // namespace powheg_dy
