#include "LEXGS.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
LEXGS<SType,Type,MeshType>::LEXGS
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    solver<SType,Type,MeshType>::smoother(dict,fvMsh),
    omega_(1.0)
{}

template<class SType, class Type, class MeshType>
void LEXGS<SType,Type,MeshType>::LEXGS::smooth
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
            meshDirection<Type,MeshType>& xd = x[d];

            const meshDirection<SType,MeshType>& Ad = A[d];
            const meshDirection<Type,MeshType>& bd = b[d];

            forAllCells(xd, i, j, k)
            {
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
