#include "zeroSigma.H"
#include "twoPhaseModel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(zeroSigma, 0);

zeroSigma::zeroSigma
(
    const twoPhaseModel& tpm,
    const dictionary& dict
)
:
    surfaceTensionScheme(tpm, dict)
{
    this->sigma_ = Zero;
}

zeroSigma::zeroSigma(const zeroSigma& sts)
:
    surfaceTensionScheme(sts)
{
    this->sigma_ = Zero;
}

zeroSigma::~zeroSigma()
{}

void zeroSigma::correct()
{
    surfaceTensionScheme::correct();
}

}

}

}
