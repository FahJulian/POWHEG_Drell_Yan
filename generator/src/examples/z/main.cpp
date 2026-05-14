#include "drell_yan.h"

using namespace powheg_dy;

int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        std::cerr << "Usage: " << argv[0] << " <config-filepath> <export-filepath>" << std::endl;
        return 1;
    }

    const std::string configPath = argv[1];
    const std::string exportPath = argv[2];

    DrellYanProcess process;

    process.init(configPath);
    process.run();
    process.writeToFile(exportPath);
}
