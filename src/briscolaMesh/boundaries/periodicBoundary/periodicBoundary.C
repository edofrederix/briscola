#include "periodicBoundary.H"
#include "addToRunTimeSelectionTable.H"
#include "level.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(periodicBoundary, 0);
addToRunTimeSelectionTable(boundary, periodicBoundary, dictionary);

periodicBoundary::periodicBoundary(const level& lvl, const dictionary& dict)
:
    parallelBoundary(lvl, dict)
{}

periodicBoundary::periodicBoundary(const periodicBoundary& pp)
:
    parallelBoundary(pp)
{}

periodicBoundary::periodicBoundary(const periodicBoundary& pp, const level& lvl)
:
    parallelBoundary(pp, lvl)
{}

}

}
