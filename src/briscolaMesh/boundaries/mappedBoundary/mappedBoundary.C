#include "mappedBoundary.H"
#include "addToRunTimeSelectionTable.H"
#include "level.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(mappedBoundary, 0);
addToRunTimeSelectionTable(boundary, mappedBoundary, dictionary);

mappedBoundary::mappedBoundary(const level& lvl, const dictionary& dict)
:
    patchBoundary(lvl, dict),
    mappingOffset_(dict.lookup("mappingOffset"))
{}

mappedBoundary::mappedBoundary(const mappedBoundary& b)
:
    patchBoundary(b),
    mappingOffset_(b.mappingOffset_)
{}

mappedBoundary::mappedBoundary(const mappedBoundary& b, const level& lvl)
:
    patchBoundary(b, lvl),
    mappingOffset_(b.mappingOffset_)
{}

}

}
