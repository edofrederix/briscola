#include "domainBoundary.H"
#include "addToRunTimeSelectionTable.H"
#include "mesh.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(domainBoundary, 0);
addToRunTimeSelectionTable(boundary, domainBoundary, dictionary);

const label domainBoundary::typeNumber = 1;

domainBoundary::domainBoundary(const mesh& msh, const dictionary& dict)
:
    boundary(msh, dict)
{}

domainBoundary::domainBoundary(const domainBoundary& pp)
:
    boundary(pp)
{}

}

}
