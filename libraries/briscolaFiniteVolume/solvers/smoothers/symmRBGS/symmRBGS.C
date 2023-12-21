#include "symmRBGS.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
symmRBGS<SType,Type,MeshType>::symmRBGS
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    solver<SType,Type,MeshType>::smoother(dict,fvMsh)
{}

template<class SType, class Type, class MeshType>
void symmRBGS<SType,Type,MeshType>::symmRBGS::smooth
(
    linearSystem<SType,Type,MeshType>& sys,
    List<Type>& xi,
    const label l,
    const label sweeps,
    const labelList& converged
)
{
    meshLevel<Type,MeshType>& x = sys.x()[l];

    const meshLevel<SType,MeshType>& A = sys.A()[l];
    const meshLevel<Type,MeshType>& b = sys.b()[l];

    bool singular = xi.size() > 0;
    xi = singular ? gAverage(x) : List<Type>(x.size(), Zero);

    for (label sweep = 0; sweep < sweeps; sweep++)
    {
        forAll(x, d)
        if (!converged[d])
        {
            meshDirection<Type,MeshType>& xd = x[d];

            const meshDirection<SType,MeshType>& Ad = A[d];
            const meshDirection<Type,MeshType>& bd = b[d];

            forAllCells(xd, i, j, k)
                if ((i+j+k)%2 == 0)
                    xd(i,j,k) =
                        (
                            bd(i,j,k)
                          - lowerRowProduct(Ad,xd,i,j,k)
                          - upperRowProduct(Ad,xd,i,j,k)
                          - xi[d]
                        )
                      / Ad(i,j,k).center();

            forAllCells(xd, i, j, k)
                if ((i+j+k)%2 == 1)
                    xd(i,j,k) =
                        (
                            bd(i,j,k)
                          - lowerRowProduct(Ad,xd,i,j,k)
                          - upperRowProduct(Ad,xd,i,j,k)
                          - xi[d]
                        )
                      / Ad(i,j,k).center();

            // Reversed direction

            forAllCellsReversed(xd, i, j, k)
                if ((i+j+k)%2 == 0)
                    xd(i,j,k) =
                        (
                            bd(i,j,k)
                          - lowerRowProduct(Ad,xd,i,j,k)
                          - upperRowProduct(Ad,xd,i,j,k)
                          - xi[d]
                        )
                      / Ad(i,j,k).center();

            forAllCellsReversed(xd, i, j, k)
                if ((i+j+k)%2 == 1)
                    xd(i,j,k) =
                        (
                            bd(i,j,k)
                          - lowerRowProduct(Ad,xd,i,j,k)
                          - upperRowProduct(Ad,xd,i,j,k)
                          - xi[d]
                        )
                      / Ad(i,j,k).center();
        }

        x.correctNonEliminatedBoundaryConditions();
    }

    x.correctEliminatedBoundaryConditions();
    x.correctUnsetBoundaryConditions();

    if (!singular)
        xi.clear();
}

}

}

}
