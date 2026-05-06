#include "file.h"

#include <sstream>
#include <fstream>

namespace powheg_dy
{
    File::File(const std::string& filePath)
        : m_filePath(filePath)
    {
    }

    std::string File::read() const
    {
        std::ifstream stream = std::ifstream(m_filePath.c_str());

        std::stringstream buffer;
        buffer << stream.rdbuf();

        buffer.seekg(0, buffer.end);
        int size = static_cast<int>(buffer.tellg());

        char* content = new char[size];

        buffer.seekg(std::ifstream::beg);
        buffer.read(content, size);

        stream.close();
        return content;
    }

    void File::write(const std::string& text)
    {
        std::ofstream stream = std::ofstream(m_filePath.c_str());

        stream << text;

        stream.close();
    }
    
} // namespace powheg_dy
