#include "process.h"

#include "powheg_dy/math/rand.h"
#include "powheg_dy/util/file.h"
#include "powheg_dy/emission/emission.h"
#include "powheg_dy/born_event/born_event_generator.h"
#include "powheg_dy/les_houches/les_houches_serializer.h"

#include <iostream>

namespace powheg_dy
{
    namespace 
    {
        static constexpr int __N_ACCEPTED_EVENTS = 20e3;
        static constexpr int __N_TRIAL_EVENTS = std::min(2e4, 1.5 * __N_ACCEPTED_EVENTS);
        static constexpr double __SECURITY_FACTOR = 1.1;

    } // namespace
    
    void Process::init(const std::string& pdfDataLocation, const std::string& pdfSet)
    {
        LHAPDF::setPaths(pdfDataLocation);
        m_pdfs = std::unique_ptr<LHAPDF::PDF>(LHAPDF::mkPDF(pdfSet, 0));
        
        m_bornPhSp = std::make_unique<BornPhaseSpace>(*this);
        m_realPhSp = std::make_unique<FKSRealPhaseSpace>(*this);
        m_emissionGenerator = std::make_unique<EmissionGenerator>(*this);
        m_bornGenerator = std::make_unique<BornEventGenerator>(*this);
    }

    void Process::run()
    {
        _clear();

        _determineMaxWeight();
        
        while (m_events.size() < __N_ACCEPTED_EVENTS)
        {   
            m_nEventTrials++;
            
            double rands[3] = { rand(), rand(), rand() };
            BornPhSpPt born = m_bornPhSp->samplePoint(rands);
            m_bornPhSp->reconstructMomenta(born);
            m_bornGenerator->computeWeightAndSampleChannel(born);
            
            assert(born.weight <= m_maxWeight);
            
            // Unweight: Accept the event at the rate of it's weight over the reference weight
            double u = rand(0.0, 1.0);
            if (u < born.weight / m_maxWeight)
            {   
                if (bornOnly())
                {
                    Emission emission = Emission().reject();
                    RealPhSpPt real = m_realPhSp->reconstruct(born, emission.rad);

                    m_events.push_back({ born, real, emission });
                }
                else
                {
                    Emission highestPtEm = Emission().reject();

                    // Highest bid procedure, for Drell Yan only ISR radiation 
                    for (int region = 1; region < 2; region++)
                    {
                        Emission emission = m_emissionGenerator->generateEmission(born, region);
                        if (emission.kt2 > highestPtEm.kt2)
                            highestPtEm = emission;
                    }
                    
                    RealPhSpPt real = m_realPhSp->reconstruct(born, highestPtEm.rad);
                    m_events.push_back({ born, real, highestPtEm });
                }

                if (m_events.size() % 1000 == 0)
                    std::cout << m_events.size() << " Events generated." << std::endl;
            }
        }

        double rejected = 0.0;
        for (const auto& event : m_events)
        {
            if (event.emission.rejected)
                rejected++;
        }

        std::cout << "No emission probability: " << rejected / m_events.size() << std::endl;
                
        std::cout << "Acceptance ratio: " << double(__N_ACCEPTED_EVENTS) / m_nEventTrials << std::endl;
                
        _computeTotalCrossSection();
    
        std::cout << "Total cross section: " << m_totalCrossSection << " pb." << std::endl;
    }

    std::string toString(const Event& event)
    {
        FourVector sum = event.real.pLMinus + event.real.pLPlus;
        double pT = sqrt(sum.pX*sum.pX + sum.pY*sum.pY);

        return std::to_string(event.born.mB) + ", "
            + std::to_string(event.born.cosTh) + ", "
            + std::to_string(pT);
    }

    void Process::writeToFile(const std::string& filePath) const
    {
        std::string fileContent;
            
        for (const auto& event : m_events)
            fileContent.append(toString(event) + '\n');
        
        File file = File(filePath);
        file.write(fileContent);

        // LesHouchesSerializer(*this).serialize(filePath);
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
            BornPhSpPt point = m_bornPhSp->samplePoint(rands);
            m_bornPhSp->reconstructMomenta(point);
            m_bornGenerator->computeWeightAndSampleChannel(point);

            if (point.weight > max_dSigma)
                max_dSigma = point.weight;
        }

        m_maxWeight = __SECURITY_FACTOR * max_dSigma;
    }

    void Process::_computeTotalCrossSection()
    {
        // const double e2 = 4.0 * PI * ALPHA();
        // * e2 * e2 / (4.0 * process.NC());
        m_totalCrossSection = m_maxWeight * __N_ACCEPTED_EVENTS / m_nEventTrials * GEV2_TO_PB();
    }

} // namespace powheg_dy
