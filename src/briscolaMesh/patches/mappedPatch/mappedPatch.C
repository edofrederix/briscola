#include "mappedPatch.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(mappedPatch, 0);
addToRunTimeSelectionTable(patch, mappedPatch, dictionary);

mappedPatch::mappedPatch
(
    const geometry& g,
    const label num,
    const word name,
    const dictionary& dict
)
:
    patch(g, num, name, dict)
{}

mappedPatch::mappedPatch(const mappedPatch& p)
:
    patch(p)
{}

mappedPatch::~mappedPatch()
{}

}

}
