#include "emptyPartPatch.H"
#include "addToRunTimeSelectionTable.H"
#include "mesh.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(emptyPartPatch, 0);
addToRunTimeSelectionTable(partPatch, emptyPartPatch, dictionary);

const label emptyPartPatch::typeNumber = 0;

emptyPartPatch::emptyPartPatch(const mesh& msh, const dictionary& dict)
:
    partPatch(msh, dict)
{
    // Empty patches are extended patches

    this->extend();
}

emptyPartPatch::emptyPartPatch(const emptyPartPatch& pp)
:
    partPatch(pp)
{}

}

}
