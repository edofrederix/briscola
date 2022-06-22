#include "parallelPartPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "mesh.H"

namespace Foam
{

namespace briscola
{

namespace partPatches
{

defineTypeNameAndDebug(parallelPartPatch, 0);
addToRunTimeSelectionTable(partPatch, parallelPartPatch, dictionary);

parallelPartPatch::parallelPartPatch(const mesh& msh, const dictionary& dict)
:
    partPatch(msh, dict),
    neighborProcNum_(readLabel(dict.lookup("neighborProcNum")))
{}

}

}

}
