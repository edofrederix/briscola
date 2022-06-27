#include "periodicPartPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "mesh.H"

namespace Foam
{

namespace briscola
{

namespace partPatches
{

defineTypeNameAndDebug(periodicPartPatch, 0);
addToRunTimeSelectionTable(partPatch, periodicPartPatch, dictionary);

periodicPartPatch::periodicPartPatch(const mesh& msh, const dictionary& dict)
:
    parallelPartPatch(msh, dict)
{}

}

}

}
