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
{
    const labelVector bo(dict.lookup("boundaryOffset"));
    const labelVector no(dict.lookup("neighborOffset"));

    const label n1 =
        cmptSum(cmptMag(bo)) == 1
      ? faceNumber(bo)
      : cmptSum(cmptMag(bo)) == 2
      ? edgeNumber(bo)
      : vertexNumber(bo);

    const label n2 =
        cmptSum(cmptMag(bo)) == 1
      ? faceNumber(bo)
      : cmptSum(cmptMag(bo)) == 2
      ? edgeNumber(bo)
      : vertexNumber(bo);

    if (n1 % 2 == 0 && n2 % 2 == 1)
    {
        master_ = true;
    }
    else if (n1 % 2 == 1 && n2 % 2 == 0)
    {
        master_ = false;
    }
    else
    {
        master_ = Pstream::myProcNo() < neighborProcNum_;
    }
}

}

}

}
