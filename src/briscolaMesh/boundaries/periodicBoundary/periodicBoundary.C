#include "periodicBoundary.H"
#include "addToRunTimeSelectionTable.H"
#include "mesh.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(periodicBoundary, 0);
addToRunTimeSelectionTable(boundary, periodicBoundary, dictionary);

const label periodicBoundary::typeNumber = 4;

periodicBoundary::periodicBoundary(const mesh& msh, const dictionary& dict)
:
    parallelBoundary(msh, dict)
{}

periodicBoundary::periodicBoundary(const periodicBoundary& pp)
:
    parallelBoundary(pp)
{}

}

}
