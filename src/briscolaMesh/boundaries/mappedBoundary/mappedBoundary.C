#include "mappedBoundary.H"
#include "addToRunTimeSelectionTable.H"
#include "mesh.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(mappedBoundary, 0);
addToRunTimeSelectionTable(boundary, mappedBoundary, dictionary);

const label mappedBoundary::typeNumber = 2;

mappedBoundary::mappedBoundary(const mesh& msh, const dictionary& dict)
:
    domainBoundary(msh, dict),
    mappingOffset_(dict.lookup("mappingOffset"))
{}

mappedBoundary::mappedBoundary(const mappedBoundary& b)
:
    domainBoundary(b),
    mappingOffset_(b.mappingOffset_)
{}

}

}
