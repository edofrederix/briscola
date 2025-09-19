#include "mappedBoundary.H"
#include "addToRunTimeSelectionTable.H"
#include "mesh.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(mappedBoundary, 0);
addToRunTimeSelectionTable(boundary, mappedBoundary, dictionary);

const label mappedBoundary::typeNumber = 1;

mappedBoundary::mappedBoundary(const mesh& msh, const dictionary& dict)
:
    patchBoundary(msh, dict),
    mappingOffset_(dict.lookup("mappingOffset"))
{}

mappedBoundary::mappedBoundary(const mappedBoundary& b)
:
    patchBoundary(b),
    mappingOffset_(b.mappingOffset_)
{}

}

}
