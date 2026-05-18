#include "les_houches_serializer.h"

#include "powheg/event.h"
#include "powheg/process.h"

namespace powheg
{
namespace 
{
    std::array<char, 8 * sizeof(long long)> STRING_CAST_BUFFER;

    std::string doubleToString(double value, int precision, bool scientific)
    {
        auto result = std::to_chars(STRING_CAST_BUFFER.begin(), STRING_CAST_BUFFER.end(),
            value, scientific ? std::chars_format::scientific : std::chars_format::fixed, precision);
            
        std::string str = std::string(STRING_CAST_BUFFER.begin(), result.ptr - STRING_CAST_BUFFER.begin());
        return value < 0 ? str : " " + str;
    }

    void writeHeader(std::ostream& content)
    {
        content <<
            "<LesHouchesEvents version=\"3.0\">\n"
            "<header>\n"
            "    <generator name=\"powheg\">Toy DY generator</generator>\n"
            "</header>\n";
    }

    void writeInitBlock(std::ostream& content, double eBeam, double sigma)
    {
        content
            << "<init>\n"
            << "     2212     2212  " 
            << doubleToString(eBeam, 5, true) << "  " 
            << doubleToString(eBeam, 5, true) 
            << "     -1     -1     -1     -1     -4      1\n"
            << "    " << doubleToString(sigma, 5, 0) << "     "
            << doubleToString(0.0, 5, true) << "  " 
            << doubleToString(1.0, 5, true) << "  " 
            << "10011\n"
            << "</init>\n";
    }

    void writeParticle(std::ostream& content, int id, int status, int mother1, 
        int mother2, int color1, int color2, const FourVector& p, double mass, double spin = 9.0)
    {
        std::string idString = (id < -9) ? std::to_string(id) 
            : (id < 0) ? " " + std::to_string(id) 
            : (id < 10) ? "  " + std::to_string(id) 
            : " " + std::to_string(id);
        std::string statusString = status < 0 ? std::to_string(status) : " " + std::to_string(status);
        std::string color1String = color1 == 0 ? "  0" : std::to_string(color1);
        std::string color2String = color2 == 0 ? "  0" : std::to_string(color2);

        content 
            << "     " << idString << "    " << statusString << "     " 
            << mother1 << "     " << mother2 << "   " 
            << color1String << "   " << color2String << "  "
            << doubleToString(p.x, 9, true) << "  " << doubleToString(p.y, 9, true) << "  " 
            << doubleToString(p.z, 9, true) << "  " << doubleToString(p.e, 9, true) << "  "
            << doubleToString(mass, 9, true) << " " << doubleToString(0.0, 5, true) << "  "
            << doubleToString(spin, 3, true) << "\n";
    }

} // anonymous namespace

    void LesHouchesSerializer::serialize(const std::string& filePath)
    {
        std::ofstream file(filePath);

        if (!file)
        {
            Log::err << "Can't write LHE file: Could not open file " << filePath << Log::endl;
            return;
        }
        
        writeHeader(file);
        writeInitBlock(file, m_config.SQRT_S / 2.0, m_process.getSigma());

        for (const auto& event : m_process.getEvents())
            writeEvent(event, file);
        
        file << "</LesHouchesEvents>\n";

    }

    void LesHouchesSerializer::writeEventHeader(std::ostream& content, int sign, int nParticles, double scalup) const
    {
        // Number of particles, process label, weight of the event
        content << "<event>\n" 
            << "       " << nParticles << " 10011  " 
            << doubleToString(static_cast<double>(sign), 5, true) << "  "
            << doubleToString(scalup, 5, true) << "  "
            // << doubleToString(m_config.ALPHA_EW, 5, true) << "  " 
            << doubleToString(-1.0, 5, true) << "  " 
            << doubleToString(m_config.alphaSCMW(scalup * scalup), 5, true) << "\n";
            // << doubleToString(m_config.alphaS(alphaSScale * alphaSScale), 5, true) << "\n";
    }

    void LesHouchesSerializer::writeEvent(const Event& event, std::ostream& content) const
    {
        if (event.emission.rejected)
        {
            writeEventBorn(event, content);
        }
        else 
        {
            if (event.emission.channel.outIDs[2] == 21)
                writeEventqqbar(event, content);
            else if (event.emission.channel.id1 == 21)
                writeEventGluonLeg1(event, content);
            else 
                writeEventGluonLeg2(event, content);
        }
    }

    void LesHouchesSerializer::writeEventBorn(const Event& event, std::ostream& content) const
    {
        const int color = 501;
        const int anticolor = 501;

        // No hardest event below the cutoff was generated -> Set SCALUP to the cutoff
        writeEventHeader(content, event.weightSign, 5, std::sqrt(m_config.PT_SQ_CUTOFF));

        if (event.bornChannel.id1 > 0)     // quark on leg 1
        {
            writeParticle(content, event.bornChannel.id1, -1, 0, 0, color, 0, event.born.p1Bar, 0.0);
            writeParticle(content, event.bornChannel.id2, -1, 0, 0, 0, anticolor, event.born.p2Bar, 0.0);
        }
        else                                // antiquark on leg 1
        {
            writeParticle(content, event.bornChannel.id1, -1, 0, 0, 0, anticolor, event.born.p1Bar, 0.0);
            writeParticle(content, event.bornChannel.id2, -1, 0, 0, color, 0, event.born.p2Bar, 0.0);
        }

        // writeParticle(content, 23, 2, 1, 2, 0, 0, event.real.pBoson, event.real.pBoson.square());

        writeParticle(content, 11,  1, 3, 3, 0, 0, event.real.pOut[0], 0.0);
        writeParticle(content, -11, 1, 3, 3, 0, 0, event.real.pOut[1], 0.0);

        content << "</event>\n";
    }

    void LesHouchesSerializer::writeEventqqbar(const Event& event, std::ostream& content) const
    {
        const int color1 = 501;
        const int color2 = 511;

        writeEventHeader(content, event.weightSign, 6, std::sqrt(event.emission.kt2));

        if (event.emission.channel.id1 > 0)     // quark on leg 1
        {
            writeParticle(content, event.emission.channel.id1, -1, 0, 0, color1, 0, event.real.p1In, 0.0);
            writeParticle(content, event.emission.channel.id2, -1, 0, 0, 0, color2, event.real.p2In, 0.0);
        }
        else                                    // antiquark on leg 1
        {
            writeParticle(content, event.emission.channel.id1, -1, 0, 0, 0, color2, event.real.p1In, 0.0);
            writeParticle(content, event.emission.channel.id2, -1, 0, 0, color1, 0, event.real.p2In, 0.0);
        }

        // writeParticle(content, 23, 2, 1, 2, 0, 0, event.real.pBoson, event.real.pBoson.square());

        writeParticle(content, 11,  1, 3, 3, 0, 0, event.real.pOut[0], 0.0);
        writeParticle(content, -11, 1, 3, 3, 0, 0, event.real.pOut[1], 0.0);

        writeParticle(content, event.emission.channel.outIDs[2], 1, 1, 2, color1, color2, event.real.pOut[2], 0.0);

        content << "</event>\n";
    }

    void LesHouchesSerializer::writeEventGluonLeg1(const Event& event, std::ostream& content) const
    {
        const int colorIn = 501;
        const int colorOut = 511;
        
        writeEventHeader(content, event.weightSign, 6, std::sqrt(event.emission.kt2));

        if (event.emission.channel.id2 > 0)     // quark on leg 2
        {
            writeParticle(content, event.emission.channel.id1, -1, 0, 0, colorOut, colorIn, event.real.p1In, 0.0);
            writeParticle(content, event.emission.channel.id2, -1, 0, 0, colorIn, 0, event.real.p2In, 0.0);
        }
        else                                    // antiquark on leg 2
        {
            writeParticle(content, event.emission.channel.id1, -1, 0, 0, colorIn, colorOut, event.real.p1In, 0.0);
            writeParticle(content, event.emission.channel.id2, -1, 0, 0, 0, colorIn, event.real.p2In, 0.0);
        }

        // writeParticle(content, 23, 2, 1, 2, 0, 0, event.real.pBoson, event.real.pBoson.square());

        writeParticle(content, 11,  1, 3, 3, 0, 0, event.real.pOut[0], 0.0);
        writeParticle(content, -11, 1, 3, 3, 0, 0, event.real.pOut[1], 0.0);

        if (event.emission.channel.id2 > 0)     // Final state quark
            writeParticle(content, event.emission.channel.outIDs[2], 1, 1, 2, colorOut, 0, event.real.pOut[2], 0.0);
        else                                    // Final state antiquark
            writeParticle(content, event.emission.channel.outIDs[2], 1, 1, 2, 0, colorOut, event.real.pOut[2], 0.0);

        content << "</event>\n";
    }

    void LesHouchesSerializer::writeEventGluonLeg2(const Event& event, std::ostream& content) const
    {
        const int colorIn = 501;
        const int colorOut = 511;
        
        writeEventHeader(content, event.weightSign, 6, std::sqrt(event.emission.kt2));

        if (event.emission.channel.id1 > 0)     // quark on leg 1
        {
            writeParticle(content, event.emission.channel.id1, -1, 0, 0, colorIn, 0, event.real.p1In, 0.0);
            writeParticle(content, event.emission.channel.id2, -1, 0, 0, colorOut, colorIn, event.real.p2In, 0.0);
        }
        else                                    // antiquark on leg 1
        {
            writeParticle(content, event.emission.channel.id1, -1, 0, 0, 0, colorIn, event.real.p1In, 0.0);
            writeParticle(content, event.emission.channel.id2, -1, 0, 0, colorIn, colorOut, event.real.p2In, 0.0);
        }

        // writeParticle(content, 23, 2, 1, 2, 0, 0, event.real.pBoson, event.real.pBoson.square());

        writeParticle(content, 11,  1, 3, 3, 0, 0, event.real.pOut[0], 0.0);
        writeParticle(content, -11, 1, 3, 3, 0, 0, event.real.pOut[1], 0.0);

        if (event.emission.channel.id1 > 0)     // Final state quark
            writeParticle(content, event.emission.channel.outIDs[2], 1, 1, 2, colorOut, 0, event.real.pOut[2], 0.0);
        else                                    // Final state antiquark
            writeParticle(content, event.emission.channel.outIDs[2], 1, 1, 2, 0, colorOut, event.real.pOut[2], 0.0);

        content << "</event>\n";
    }
    
} // namespace powheg
