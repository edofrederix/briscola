#include "uniformMesh.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(uniformMesh, 0);
addToRunTimeSelectionTable(mesh, uniformMesh, dictionary);

void uniformMesh::setCellSize()
{
    // Take the first value of the rectilinear cell size arrays

    cellSize_ = vector
    (
        cellSizes()[0][0],
        cellSizes()[1][0],
        cellSizes()[2][0]
    );
}

uniformMesh::uniformMesh(const IOdictionary& dict)
:
    rectilinearMesh(dict)
{
    setCellSize();
}

uniformMesh::uniformMesh(autoPtr<mesh>& mshPtr)
:
    rectilinearMesh(mshPtr)
{
    setCellSize();
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
