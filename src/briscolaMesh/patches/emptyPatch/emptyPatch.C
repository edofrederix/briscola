#include "emptyPatch.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(emptyPatch, 0);
addToRunTimeSelectionTable(patch, emptyPatch, dictionary);

emptyPatch::emptyPatch
(
    const geometry& g,
    const label num,
    const word name,
    const dictionary& dict
)
:
    patch(g, num, name, dict)
{}

emptyPatch::emptyPatch(const emptyPatch& p)
:
    patch(p)
{}

emptyPatch::~emptyPatch()
{}

}

}
