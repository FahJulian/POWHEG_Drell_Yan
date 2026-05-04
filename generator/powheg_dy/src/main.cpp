#include "./powheg_dy/process.h"
#include "./powheg_dy/math/rand.h"
#include "./powheg_dy/phase_space/real_phase_space.h"

const std::string exportFilePath = "/home/julian/documents/uni/master/master_thesis/learning/powheg_drell_yan/data/events/events.lhe";
const std::string pdfDataLocation = "/home/julian/documents/uni/master/master_thesis/learning/powheg_drell_yan/data/lhapdf";
const std::string pdfSet = "NNPDF40_lo_as_01180";

class DYProcess : public powheg_dy::Process
{   
public:
    inline double pt2Cutoff() const override { return 5.0; }     // GeV^2
    inline int nBornLegs() const override { return 4; }
    inline int nRealLegs() const override { return 5; }
    inline bool bornOnly() const override { return false; }
    inline double mMin() const override { return 75; }
    inline double mMax() const override { return 105; }
    inline double sqrtS() const override { return 8.0e3; }
    inline double S() const override { return 8.0e3*8.0e3; }
};

int main()
{
    DYProcess process;

    process.init(pdfDataLocation, pdfSet);
    process.run();
    process.writeToFile(exportFilePath);
}
