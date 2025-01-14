#include "parallelBoundary.H"
#include "addToRunTimeSelectionTable.H"
#include "mesh.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(parallelBoundary, 0);
addToRunTimeSelectionTable(boundary, parallelBoundary, dictionary);

const label parallelBoundary::typeNumber = 2;

parallelBoundary::parallelBoundary(const mesh& msh, const dictionary& dict)
:
    boundary(msh, dict),
    neighborProcNum_(readLabel(dict.lookup("neighborProcNum"))),
    tag_(0)
{
    // Set master bool. This boundary is master if its face, edge or vertex
    // number is even and the neighboring one is odd. If they are both odd or
    // both even, the lower processor number is master.

    const labelVector bo(dict.lookup("offset"));
    const labelVector no(dict.lookup("neighborOffset"));

    const label n1 =
        cmptSum(cmptMag(bo)) == 1
      ? faceNumber(bo)
      : cmptSum(cmptMag(bo)) == 2
      ? edgeNumber(bo)
      : vertexNumber(bo);

    const label n2 =
        cmptSum(cmptMag(no)) == 1
      ? faceNumber(no)
      : cmptSum(cmptMag(no)) == 2
      ? edgeNumber(no)
      : vertexNumber(no);

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

    // Parallel patches are extended patches

    this->extend();
}

parallelBoundary::parallelBoundary(const parallelBoundary& pp)
:
    boundary(pp),
    neighborProcNum_(pp.neighborProcNum_)
{}

}

}
