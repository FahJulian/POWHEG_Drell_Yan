#include "process.h"

#include "powheg_dy/math/rand.h"
#include "powheg_dy/util/file.h"
#include "powheg_dy/emission/emission.h"
#include "powheg_dy/born_event/born_event_generator.h"
#include "powheg_dy/les_houches/les_houches_serializer.h"

#include <iostream>
#include <sstream>

namespace powheg_dy
{
namespace 
{
    static constexpr double GEV_SQ_TO_PB = 0.389379338e9;
    static constexpr double SECURITY_FACTOR = 1.1;

} // anonymous namespace

    Process::~Process()
    {
        LHAPDF::setVerbosity(0);
        Log::info("Using LHAPDF. Please make sure to cite the paper:");
        Log::info("Eur.Phys.J. C75 (2015) 3, 132  (http://arxiv.org/abs/1412.7420)");
    }
    
    void Process::init(const std::string& pdfSet)
    {
        Log::info("Starting intitialization");

        std::stringstream buffer;

        auto* oldCout = std::cout.rdbuf(buffer.rdbuf());
        auto* oldCerr = std::cerr.rdbuf(buffer.rdbuf());

        m_config.PDF = std::unique_ptr<LHAPDF::PDF>(LHAPDF::mkPDF(pdfSet, 0));

        std::cout.rdbuf(oldCout);
        std::cerr.rdbuf(oldCerr);

        std::string line;
        while (std::getline(buffer, line, '\n'))
            Log::info << line << std::endl;

        setDependentParams(m_config);
        extractLambdaFromPdf(m_config);

        m_bornPhSp = std::make_shared<BornPhaseSpace>(m_config);
        m_realPhSp = std::make_shared<FKSRealPhaseSpace>(m_config);
        m_bornGenerator = std::make_shared<BornEventGenerator>(*this, m_config, m_bornPhSp);
        m_emissionGenerator = std::make_shared<EmissionGenerator>(*this, m_config, m_realPhSp);

        Log::info << "Intitialization complete" << std::endl << std::endl;
    }

    void Process::run()
    {
        try
        {
            clear();
            determineMaxWeight();
            generateEvents();
            computeTotalCrossSection();
            analyse();
        }
        catch(const std::runtime_error& e)
        {
            Log::err << "Aborting..." << std::endl;
        }
    }

    std::string toString(const Event& event)
    {
        FourVector sum = event.real.pLMinus + event.real.pLPlus;
        double pT = sqrt(sum.pX*sum.pX + sum.pY*sum.pY);

        return std::to_string(event.born.mBoson) + ", "
            + std::to_string(event.born.cosTh) + ", "
            + std::to_string(pT);
    }

    void Process::writeToFile(const std::string& filePath) const
    {
        // std::string fileContent;
            
        // for (const auto& event : m_events)
        //     fileContent.append(toString(event) + '\n');
        
        // File file = File(filePath);
        // file.write(fileContent);

        Log::info << "Generating LHE file at " << filePath << std::endl;
        LesHouchesSerializer(*this, m_config).serialize(filePath);
        Log::info("Done generating LHE file");
    }

    void Process::clear()
    {
        m_events.clear();
        m_events.reserve(m_config.N_ACCEPTED_EVENTS);
        m_nEventTrials = 0;
        m_totalCrossSection = 0.0;
        m_maxWeight = 0.0;
    }

    void Process::determineMaxWeight()
    {
        Log::info("Determining Born Veto weight");

        double max_dSigma = 0.0;
        
        for (int i = 0; i < m_config.N_TRIAL_EVENTS; i++)
        {   
            double rands[3] = { rand(), rand(), rand() };
            BornPhSpPt point = m_bornPhSp->samplePoint(rands);
            m_bornGenerator->computeWeightAndSampleChannel(point);

            if (point.weight > max_dSigma)
                max_dSigma = point.weight;
        }

        m_maxWeight = SECURITY_FACTOR * max_dSigma;

        Log::info << "Done determining Born Veto weight." << std::endl << std::endl;
    }

    void Process::generateEvents()  
    {
        Log::info("Starting event generation");
        while (static_cast<int>(m_events.size()) < m_config.N_ACCEPTED_EVENTS)
        {   
            m_nEventTrials++;
            
            double rands[3] = { rand(), rand(), rand() };
            BornPhSpPt born = m_bornPhSp->samplePoint(rands);
            m_bornGenerator->computeWeightAndSampleChannel(born);
            
            assert(born.weight <= m_maxWeight, "Born weight " << born.weight << " exceeds max weight " << m_maxWeight);
            
            // Unweight: Accept the event at the rate of it's weight over the reference weight
            double u = rand(0.0, 1.0);
            if (u < born.weight / m_maxWeight)
            {   
                if (m_config.BORNONLY)
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
                    Log::info << m_events.size() << " Events generated." << std::endl;
            }
        }

        Log::info << "Event generation done." << std::endl << std::endl;
    }

    void Process::computeTotalCrossSection()
    {
        m_totalCrossSection = m_maxWeight * GEV_SQ_TO_PB * 
            m_config.N_ACCEPTED_EVENTS / m_nEventTrials;
    }

    void Process::analyse()
    {
        double rejected = 0.0;
        for (const auto& event : m_events)
        {
            if (event.emission.rejected)
                rejected++;
        }

        Log::info << "No emission probability: " << rejected / static_cast<double>(m_events.size()) << std::endl;
        Log::info << "Acceptance ratio: " << double(m_config.N_ACCEPTED_EVENTS) / m_nEventTrials << std::endl;
        Log::info << "Total cross section: " << m_totalCrossSection << " pb." << std::endl;
    }

} // namespace powheg_dy
