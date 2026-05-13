#pragma once
 
#include "powheg_dy/base.h"

namespace powheg_dy
{
    class File
    {
    public:
        File(const std::string& filePath);

        ~File();

        std::string read() const;

        // std::vector<std::string> readLines() const;

        bool getLine(std::string& line);
        
        void write(const std::string& text);

    private:
        std::string m_filePath;
        std::ifstream m_stream;
    };

} // namespace powheg_dy
