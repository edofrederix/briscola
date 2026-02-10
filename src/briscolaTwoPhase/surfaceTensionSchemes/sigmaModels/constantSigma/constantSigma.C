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
    const fvMesh& fvMsh,
    const dictionary& dict,
    const normalScheme& normal,
    const colocatedScalarField& alpha
)
:
    surfaceTensionScheme(fvMsh, dict, normal, alpha)
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
