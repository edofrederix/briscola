#include "uniformMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "parallelBoundary.H"

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

    const PtrList<PartialList<scalar>>& x = localPoints();

    // Check if the point is outside the domain boundaries

    for (label d = 0; d < 3; d++)
        if (q[d] < x[d].first() - 1e-12 || q[d] > x[d].last() + 1e-12)
            return -unitXYZ;

    // Check if the point is just above/below the upper/lower boundary. If so,
    // the point is included if the boundary is not be a parallel one. If the
    // point is exactly on the lower boundary it must be included anyway.

    for (label d = 0; d < 3; d++)
        if (q[d] >= x[d].last() && q[d] <= x[d].last() + 1e-12)
            if (b(units[d]).castable<parallelBoundary>())
                return -unitXYZ;

    for (label d = 0; d < 3; d++)
        if (q[d] < x[d].first() && q[d] >= x[d].first() - 1e-12)
            if (b(-units[d]).castable<parallelBoundary>())
                return -unitXYZ;

    // Algebraic relation. Cast to label automatically applies floor. Assign
    // points that are on the upper boundary to the nearest internal cell.

    const labelVector N(this->operator[](0).N());
    const labelVector R
    (
        cmptDivide(N, this->operator[](l).N())
    );

    labelVector ijk;

    for (label d = 0; d < 3; d++)
        ijk[d] =
            Foam::min
            (
                Foam::max
                (
                    (q[d] - x[d].first())/(x[d].last() - x[d].first())*N[d],
                    0
                ),
                N[d]-1
            );

    return cmptDivide(ijk, R);
}


}

}
