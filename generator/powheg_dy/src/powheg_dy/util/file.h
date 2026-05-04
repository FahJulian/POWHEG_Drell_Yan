#pragma once
 

#include <string>

namespace powheg_dy
{
    class File
    {
    public:
        File(const std::string& filePath);

        std::string read() const;

        // std::vector<std::string> readLines() const;

        void write(const std::string& text);

    private:
        std::string m_filePath;
    };

} // namespace powheg_dy
