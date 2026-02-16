#include "patchBoundary.H"
#include "addToRunTimeSelectionTable.H"
#include "level.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(patchBoundary, 0);
addToRunTimeSelectionTable(boundary, patchBoundary, dictionary);

patchBoundary::patchBoundary(const level& lvl, const dictionary& dict)
:
    boundary(lvl, dict)
{}

patchBoundary::patchBoundary(const patchBoundary& pp)
:
    boundary(pp)
{}

patchBoundary::patchBoundary(const patchBoundary& pp, const level& lvl)
:
    boundary(pp, lvl)
{}

}

}
