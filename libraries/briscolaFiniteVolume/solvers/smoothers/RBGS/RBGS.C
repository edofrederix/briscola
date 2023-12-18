#include "RBGS.H"
#include "directSolver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
RBGS<SType,Type,MeshType>::RBGS
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    solver<SType,Type,MeshType>::smoother(dict,fvMsh),
    omega_(1.0)
{}

template<class SType, class Type, class MeshType>
void RBGS<SType,Type,MeshType>::RBGS::smooth
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

    if (!singular)
        xi = List<Type>(x.size(), Zero);

    for (label sweep = 0; sweep < sweeps; sweep++)
    {
        if (singular)
            forAll(xi, d)
                xi[d] = gAverage(x[d]);

        forAll(x, d)
        if (!converged[d])
        {
            // Even points

            meshDirection<Type,MeshType>& xd = x[d];

            const meshDirection<SType,MeshType>& Ad = A[d];
            const meshDirection<Type,MeshType>& bd = b[d];

            forAllCells(xd, i, j, k)
            {
                if ((i+j+k) % 2 == 0)
                    xd(i,j,k) +=
                        omega_
                      * (bd(i,j,k) - rowProduct(Ad,xd,i,j,k) - xi[d])
                      / Ad(i,j,k).center();
            }
        }

        x.correctCommBoundaryConditions();

        forAll(x, d)
        if (!converged[d])
        {
            // Odd points

            meshDirection<Type,MeshType>& xd = x[d];

            const meshDirection<SType,MeshType>& Ad = A[d];
            const meshDirection<Type,MeshType>& bd = b[d];

            forAllCells(xd, i, j, k)
            {
                if ((i+j+k) % 2 == 1)
                    xd(i,j,k) +=
                        omega_
                      * (bd(i,j,k) - rowProduct(Ad,xd,i,j,k) - xi[d])
                      / Ad(i,j,k).center();
            }
        }

        x.correctCommBoundaryConditions();
    }

    x.correctNonCommBoundaryConditions();

    if (!singular)
        xi.clear();
}

}

}

}
