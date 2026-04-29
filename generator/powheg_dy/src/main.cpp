#include "./powheg_dy/process.h"

constexpr double sqrtS = 8.0e3; // GeV
constexpr double mMin = 20; // GeV
constexpr double mMax = 200; // GeV
const std::string exportFilePath = "/home/julian/documents/uni/master/master_thesis/learning/simple_drell_yan/data/events/events.lhe";
const std::string pdfDataLocation = "/home/julian/documents/uni/master/master_thesis/learning/simple_drell_yan/data/lhapdf";
const std::string pdfSet = "NNPDF40_lo_as_01180";

int main()
{
    powheg_dy::Process process(sqrtS, mMin, mMax);

    process.init(pdfDataLocation, pdfSet);
    process.run();
    process.writeToFile(exportFilePath);
}
