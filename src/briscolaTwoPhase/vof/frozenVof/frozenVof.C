#include "frozenVof.H"

#include "addToRunTimeSelectionTable.H"
#include "rectilinearMesh.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(frozenVof, 0);
addToRunTimeSelectionTable(vof, frozenVof, dictionary);

frozenVof::frozenVof
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    normalScheme& normal,
    vofField& alpha
)
:
    vof(fvMsh, dict, normal, alpha)
{}

frozenVof::frozenVof(const frozenVof& vf)
:
    vof(vf)
{}

void frozenVof::solve(const colocatedScalarFaceField&)
{}

}

}

}
