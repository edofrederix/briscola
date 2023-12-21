#include "symmJAC.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
symmJAC<SType,Type,MeshType>::symmJAC
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    solver<SType,Type,MeshType>::smoother(dict,fvMsh),
    omega_(0.8)
{}

template<class SType, class Type, class MeshType>
void symmJAC<SType,Type,MeshType>::symmJAC::smooth
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

    meshLevel<Type,MeshType> y(x);

    for (label sweep = 0; sweep < sweeps; sweep++)
    {
        forAll(x, d)
        if (!converged[d])
        {
            meshDirection<Type,MeshType>& yd = y[d];

            const meshDirection<Type,MeshType>& xd = x[d];
            const meshDirection<SType,MeshType>& Ad = A[d];
            const meshDirection<Type,MeshType>& bd = b[d];

            forAllCells(xd, i, j, k)
                yd(i,j,k) =
                    yd(i,j,k)*(1.0-omega_)
                  + omega_
                  * (
                        bd(i,j,k)
                      - lowerRowProduct(Ad,xd,i,j,k)
                      - upperRowProduct(Ad,xd,i,j,k)
                      - xi[d]
                    )
                  / Ad(i,j,k).center();

            // Reversed direction

            forAllCellsReversed(xd, i, j, k)
                yd(i,j,k) =
                    yd(i,j,k)*(1.0-omega_)
                  + omega_
                  * (
                        bd(i,j,k)
                      - lowerRowProduct(Ad,yd,i,j,k)
                      - upperRowProduct(Ad,yd,i,j,k)
                      - xi[d]
                    )
                  / Ad(i,j,k).center();
        }

        x = y;

        x.correctNonEliminatedBoundaryConditions();
    }

    x.correctEliminatedBoundaryConditions();

    if (!singular)
        xi.clear();
}

}

}

}
