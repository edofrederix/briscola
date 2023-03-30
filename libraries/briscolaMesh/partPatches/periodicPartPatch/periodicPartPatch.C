#include "periodicPartPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "mesh.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(periodicPartPatch, 0);
addToRunTimeSelectionTable(partPatch, periodicPartPatch, dictionary);

const label periodicPartPatch::typeNumber = 2;

periodicPartPatch::periodicPartPatch(const mesh& msh, const dictionary& dict)
:
    parallelPartPatch(msh, dict)
{}

periodicPartPatch::periodicPartPatch(const periodicPartPatch& pp)
:
    parallelPartPatch(pp)
{}

}

}
