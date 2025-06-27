#include "geometricVof.H"

#include "addToRunTimeSelectionTable.H"
#include "rectilinearMesh.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(geometricVof, 0);

geometricVof::geometricVof
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    normalScheme& normal,
    colocatedScalarField& alpha
)
:
    vof(fvMsh, dict, normal, alpha),
    rectilinear_(fvMsh.msh()[0].rectilinear() == unitXYZ),
    lve_(rectilinear_)
{}

geometricVof::geometricVof(const geometricVof& vf)
:
    vof(vf),
    rectilinear_(vf.rectilinear_),
    lve_(rectilinear_)
{}

geometricVof::~geometricVof()
{}

}

}

}
