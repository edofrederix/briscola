#include "emptyBoundary.H"
#include "addToRunTimeSelectionTable.H"
#include "level.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(emptyBoundary, 0);
addToRunTimeSelectionTable(boundary, emptyBoundary, dictionary);

emptyBoundary::emptyBoundary(const level& lvl, const dictionary& dict)
:
    boundary(lvl, dict)
{}

emptyBoundary::emptyBoundary(const emptyBoundary& pp)
:
    boundary(pp)
{}

emptyBoundary::emptyBoundary(const emptyBoundary& pp, const level& lvl)
:
    boundary(pp, lvl)
{}

}

}
