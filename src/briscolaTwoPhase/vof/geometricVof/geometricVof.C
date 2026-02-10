#include "geometricVof.H"

#include "addToRunTimeSelectionTable.H"
#include "rectilinearMesh.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

geometricVof::geometricVof
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    normalScheme& normal,
    vofField& alpha
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

}

}

}
