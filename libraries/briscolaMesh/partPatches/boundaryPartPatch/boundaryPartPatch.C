#include "boundaryPartPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "mesh.H"

namespace Foam
{

namespace briscola
{

namespace partPatches
{

defineTypeNameAndDebug(boundaryPartPatch, 0);
addToRunTimeSelectionTable(partPatch, boundaryPartPatch, dictionary);

boundaryPartPatch::boundaryPartPatch(const mesh& msh, const dictionary& dict)
:
    partPatch(msh, dict)
{}

boundaryPartPatch::boundaryPartPatch(const boundaryPartPatch& pp)
:
    partPatch(pp)
{}

}

}

}
