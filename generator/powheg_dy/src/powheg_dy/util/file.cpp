#include "file.h"

namespace powheg_dy
{
    File::File(const std::string& filePath)
        : m_filePath(filePath), m_stream(std::ifstream(m_filePath.c_str()))
    {
    }

    File::~File()
    {
        m_stream.close();
    }

    std::string File::read() const
    {
        std::stringstream buffer;
        buffer << m_stream.rdbuf();

        buffer.seekg(0, buffer.end);
        int size = static_cast<int>(buffer.tellg());

        if (size == -1)
            return "";

        char* content = new char[size];

        buffer.seekg(std::ifstream::beg);
        buffer.read(content, size);

        return content;
    }

    bool File::getLine(std::string& line) {
        return static_cast<bool>(std::getline(m_stream, line));
    }

    void File::write(const std::string& text)
    {
        std::ofstream stream = std::ofstream(m_filePath.c_str());

        stream << text;

        stream.close();
    }
    
} // namespace powheg_dy
