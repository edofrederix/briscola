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
    const fvMesh& fvMsh,
    const dictionary& dict,
    const normalScheme& normal,
    const colocatedScalarField& alpha
)
:
    surfaceTensionScheme(fvMsh, dict, normal, alpha)
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
