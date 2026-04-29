#include "les_houches_serializer.h"

#include "powheg_dy/file.h"
#include "powheg_dy/event.h"
#include "powheg_dy/process.h"

namespace powheg_dy
{
    namespace 
    {
        void __writeHeader(std::stringstream& content)
        {
            content <<
                "<LesHouchesEvents version=\"3.0\">\n"
                "<header>\n"
                "    <generator name=\"powheg_dy\">Toy DY generator</generator>\n"
                "</header>\n";
        }

        void __writeInitBlock(std::stringstream& content, double eBeam, double sigma)
        {
            content
                << "<init>\n"
                << "    2212 2212 " << eBeam << " " << eBeam << " 0 0 0 0 3 1\n"
                << "    " << sigma << " 0.0 1.0 1001\n"
                << "</init>\n";
        }

        void __writeParticle(std::stringstream& content, int id, int status, int mother1, 
            int mother2, int color1, int color2, const FourVector& p, double mass = 0.0, double spin = 9.0)
        {
            content 
                << "    " << id << " " << status << " " 
                << mother1 << " " << mother2 << " " 
                << color1 << " " << color2 << " "
                << p.x << " " << p.y << " " << p.z << " " << p.e << " "
                << mass << " " << 0.0 << " "<< spin << "\n";
        }

        void __writeEvent(const Event& event, std::stringstream& content)
        {
            bool hasGluon = !event.getEmission().isRejected();

            double bornScale = hasGluon ? std::sqrt(event.getEmission().getT()) : event.getBornEvent().getMBoson();
            int nParticles = hasGluon ? 5 : 4;

            content << "<event>\n" 
                << "    " << nParticles << " 1001 1.0 "   // Number of particles, process label, weight of the event
                << bornScale << " "
                << Physics::ALPHA << " " << Physics::alphaSOneLoop(bornScale*bornScale, 5) << "\n";

            int color = 501;
            int anticolor = hasGluon ? 502 : 501;

            if (event.getBornEvent().getPartonId() > 0)     // quark on leg 1 
            {
                __writeParticle(content,  event.getBornEvent().getPartonId(), -1, 0, 0, color, 0, event.getP1In());
                __writeParticle(content, -event.getBornEvent().getPartonId(), -1, 0, 0, 0, anticolor, event.getP2In());
            }
            else                // antiquark on leg 1
            {
                __writeParticle(content,  event.getBornEvent().getPartonId(), -1, 0, 0, 0, anticolor, event.getP1In());
                __writeParticle(content, -event.getBornEvent().getPartonId(), -1, 0, 0, color, 0, event.getP2In());
            }

            __writeParticle(content, 13,  1, 1, 2, 0, 0, event.getP1Out());
            __writeParticle(content, -13, 1, 1, 2, 0, 0, event.getP2Out());

            if (hasGluon)
                __writeParticle(content, 21, 1, 1, 2, color, anticolor, event.getPGluon());

            content << "</event>\n";
        }

    } // namespace

    void LesHouchesSerializer::serialize(const std::string& filePath)
    {
        std::stringstream content;

        __writeHeader(content);
        __writeInitBlock(content, m_process.getSqrtS() / 2.0, m_process.getSigma());

        for (const auto& event : m_process.getEvents())
            __writeEvent(event, content);
        
        content << "</LesHouchesEvents>\n";
        
        File(filePath).write(content.str());
    }
    
} // namespace powheg_dy
