#include "noSurfaceTension.H"
#include "twoPhaseModel.H"
#include "exSchemes.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SigmaModel>
noSurfaceTension<SigmaModel>::noSurfaceTension
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    const normalScheme& normal,
    const colocatedScalarField& alpha
)
:
    SigmaModel(fvMsh, dict, normal, alpha)
{
    colocatedScalarFaceField::operator=(Zero);
}

template<class SigmaModel>
noSurfaceTension<SigmaModel>::noSurfaceTension(const noSurfaceTension& sts)
:
    SigmaModel(sts)
{
    colocatedScalarFaceField::operator=(Zero);
}

template<class SigmaModel>
noSurfaceTension<SigmaModel>::~noSurfaceTension()
{}

template<class SigmaModel>
void noSurfaceTension<SigmaModel>::correct()
{
    SigmaModel::correct();
}

}

}

}
