#include "patchBoundary.H"
#include "addToRunTimeSelectionTable.H"
#include "mesh.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(patchBoundary, 0);
addToRunTimeSelectionTable(boundary, patchBoundary, dictionary);

const label patchBoundary::typeNumber = 2;

patchBoundary::patchBoundary(const mesh& msh, const dictionary& dict)
:
    boundary(msh, dict)
{}

patchBoundary::patchBoundary(const patchBoundary& pp)
:
    boundary(pp)
{}

}

}
