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

uniformMesh::~uniformMesh()
{}

labelVector uniformMesh::findCell(const vector& p, const label l) const
{
    const vector q
    (
        trimPrecision
        (
            vector
            (
                p & base().x(),
                p & base().y(),
                p & base().z()
            )
        )
    );

    const scalarList& x = localPoints()[0];
    const scalarList& y = localPoints()[1];
    const scalarList& z = localPoints()[2];

    // If the point is on the upper boundary include it in the last cell

    if
    (
        q.x() < x[0] || q.x() > x[x.size()-1]
     || q.y() < y[0] || q.y() > y[y.size()-1]
     || q.z() < z[0] || q.z() > z[z.size()-1]
    )
    {
        return -unitXYZ;
    }

    const labelVector N(this->operator[](0).N());
    const labelVector R
    (
        cmptDivide(N,this->operator[](l).N())
    );

    // Algebraic relation. Cast to label automatically applies floor. Assign
    // points that are exactly on the upper boundary to the nearest internal
    // cell.

    const label i =
        Foam::min((q.x() - x[0])/(x[x.size()-1] - x[0])*N.x(), N.x()-1);
    const label j =
        Foam::min((q.y() - y[0])/(y[y.size()-1] - y[0])*N.y(), N.y()-1);
    const label k =
        Foam::min((q.z() - z[0])/(z[z.size()-1] - z[0])*N.z(), N.z()-1);

    return labelVector(i/R.x(), j/R.y(), k/R.z());
}


}

}
