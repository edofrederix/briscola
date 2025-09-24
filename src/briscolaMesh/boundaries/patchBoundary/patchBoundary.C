#include "patchBoundary.H"
#include "addToRunTimeSelectionTable.H"
#include "mesh.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(patchBoundary, 0);
addToRunTimeSelectionTable(boundary, patchBoundary, dictionary);

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
