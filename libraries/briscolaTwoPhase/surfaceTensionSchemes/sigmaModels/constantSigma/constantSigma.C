#include "constantSigma.H"
#include "twoPhaseModel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(constantSigma, 0);

constantSigma::constantSigma
(
    const twoPhaseModel& tpm,
    const dictionary& dict
)
:
    surfaceTensionScheme(tpm, dict)
{
    this->sigma_ = readScalar(dict.lookup("sigma"));
}

constantSigma::constantSigma(const constantSigma& sts)
:
    surfaceTensionScheme(sts)
{}

constantSigma::~constantSigma()
{}

void constantSigma::correct()
{
    surfaceTensionScheme::correct();
}

}

}

}
