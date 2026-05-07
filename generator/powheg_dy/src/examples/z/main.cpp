#include "drell_yan.h"

using namespace powheg_dy;

const std::string exportFilePath = "/home/julian/documents/uni/master/master_thesis/learning/powheg_drell_yan/data/events/events.lhe";
const std::string pdfSet = "NNPDF40MC_nlo_as_01180";

int main()
{
    DrellYanProcess process;

    process.init(pdfSet);
    process.run();
    process.writeToFile(exportFilePath);
}
