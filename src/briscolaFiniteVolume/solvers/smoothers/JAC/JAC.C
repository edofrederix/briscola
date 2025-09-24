#include "JAC.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
const scalar JAC<SType,Type,MeshType>::JAC::omega_ = 0.8;

template<class SType, class Type, class MeshType>
void JAC<SType,Type,MeshType>::JAC::Smooth
(
    linearSystem<SType,Type,MeshType>& sys,
    const label l,
    const label sweeps,
    const labelList& converged
)
{
    meshLevel<Type,MeshType>& x = sys.x()[l];

    const meshLevel<SType,MeshType>& A = sys.A()[l];
    const meshLevel<Type,MeshType>& b = sys.b()[l];
    const meshLevel<label,MeshType>& f = sys.forcingMask()[l];

    const List<bool> diagonal(sys.diagonal());

    meshLevel<Type,MeshType> y(x);

    for (label sweep = 0; sweep < sweeps; sweep++)
    {
        forAll(x, d)
        if (!converged[d])
        {
            if (diagonal[d])
            {
                solver<SType,Type,MeshType>::smoother::smoothDiag(sys, l, d);
            }
            else
            {
                meshDirection<Type,MeshType>& yd = y[d];

                const meshDirection<Type,MeshType>& xd = x[d];
                const meshDirection<SType,MeshType>& Ad = A[d];
                const meshDirection<Type,MeshType>& bd = b[d];
                const meshDirection<label,MeshType>& fd = f[d];

                forAllCells(xd, i, j, k)
                {
                    if (!fd(i,j,k))
                    {
                        yd(i,j,k) =
                            yd(i,j,k)*(1.0-omega_)
                          + omega_
                          * (
                                bd(i,j,k)
                              - offDiagRowProduct(Ad,xd,i,j,k)
                            )
                          / Ad(i,j,k).center();
                    }
                }
            }
        }

        x = y;

        if (l==0)
            x.correctImmersedBoundaryConditions();

        x.correctNonEliminatedBoundaryConditions();
    }

    x.correctEliminatedBoundaryConditions();
    x.correctUnsetBoundaryConditions();
}

template<class SType, class Type, class MeshType>
JAC<SType,Type,MeshType>::JAC
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    solver<SType,Type,MeshType>::smoother(dict,fvMsh)
{}

template<class SType, class Type, class MeshType>
void JAC<SType,Type,MeshType>::JAC::smooth
(
    linearSystem<SType,Type,MeshType>& sys,
    const label l,
    const label sweeps,
    const labelList& converged
)
{
    meshLevel<Type,MeshType>& x = sys.x()[l];

    const meshLevel<SType,MeshType>& A = sys.A()[l];
    const meshLevel<Type,MeshType>& b = sys.b()[l];
    const meshLevel<label,MeshType>& f = sys.forcingMask()[l];

    const List<bool> diagonal(sys.diagonal());

    meshLevel<Type,MeshType> y(x);

    for (label sweep = 0; sweep < sweeps; sweep++)
    {
        forAll(x, d)
        if (!converged[d])
        {
            if (diagonal[d])
            {
                this->smoothDiag(sys, l, d);
            }
            else
            {
                meshDirection<Type,MeshType>& yd = y[d];

                const meshDirection<Type,MeshType>& xd = x[d];
                const meshDirection<SType,MeshType>& Ad = A[d];
                const meshDirection<Type,MeshType>& bd = b[d];
                const meshDirection<label,MeshType>& fd = f[d];

                forAllCells(xd, i, j, k)
                {
                    if (!fd(i,j,k))
                    {
                        yd(i,j,k) =
                            yd(i,j,k)*(1.0-omega_)
                          + omega_
                          * (
                                bd(i,j,k)
                              - lowerRowProduct(Ad,xd,i,j,k)
                              - upperRowProduct(Ad,xd,i,j,k)
                            )
                          / Ad(i,j,k).center();
                    }
                }
            }
        }

        x = y;

        if (l==0)
            x.correctImmersedBoundaryConditions();

        x.correctNonEliminatedBoundaryConditions();
    }

    x.correctEliminatedBoundaryConditions();
}

}

}

}
