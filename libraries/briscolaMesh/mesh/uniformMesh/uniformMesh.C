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

labelVector uniformMesh::findCell(const vector& p, const label l) const
{
    const vector q(p & base().x(), p & base().y(), p & base().z());

    const scalarList& x = localPoints()[0];
    const scalarList& y = localPoints()[1];
    const scalarList& z = localPoints()[2];

    if
    (
        q.x() < x[0] || q.x() >= x[x.size()-1]
     || q.y() < y[0] || q.y() >= y[y.size()-1]
     || q.z() < z[0] || q.z() >= z[z.size()-1]
    )
    {
        return -unitXYZ;
    }

    const labelVector N(this->operator[](0).N());
    const labelVector R
    (
        cmptDivide(N,this->operator[](l).N())
    );

    // Algebraic relation

    const label i = (q.x() - x[0])/(x[x.size()-1] - x[0])*N.x();
    const label j = (q.y() - y[0])/(y[y.size()-1] - y[0])*N.y();
    const label k = (q.z() - z[0])/(z[z.size()-1] - z[0])*N.z();

    return labelVector(i/R.x(), j/R.y(), k/R.z());
}


}

}
