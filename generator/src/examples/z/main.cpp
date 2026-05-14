#include "drell_yan.h"

using namespace powheg_dy;

const std::string exportFilePath = "/home/julian/documents/uni/master/master_thesis/learning/powheg_drell_yan/data/events/events.lhe";
const std::string pdfSet = "NNPDF40MC_nlo_as_01180";

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

    process.init(configPath, pdfSet);
    process.run();
    process.writeToFile(exportPath);
}
