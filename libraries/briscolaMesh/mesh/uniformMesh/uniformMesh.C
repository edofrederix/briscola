#include "uniformMesh.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(uniformMesh, 0);
addToRunTimeSelectionTable(mesh, uniformMesh, dictionary);

void uniformMesh::setMetrics()
{
    // Take the first value of the rectilinear cell size arrays

    cellSize_ = vector
    (
        localCellSizes()[0][0],
        localCellSizes()[1][0],
        localCellSizes()[2][0]
    );
}

uniformMesh::uniformMesh(const IOdictionary& dict)
:
    rectilinearMesh(dict)
{
    setMetrics();
}

uniformMesh::uniformMesh(autoPtr<mesh>& mshPtr)
:
    rectilinearMesh(mshPtr)
{
    setMetrics();
}

uniformMesh::uniformMesh(const uniformMesh& msh)
:
    rectilinearMesh(msh),
    cellSize_(msh.cellSize_)
{}

uniformMesh::uniformMesh(uniformMesh& msh, bool reuse)
:
    rectilinearMesh(msh, reuse),
    cellSize_(msh.cellSize_)
{}

uniformMesh::~uniformMesh()
{}

}

}
