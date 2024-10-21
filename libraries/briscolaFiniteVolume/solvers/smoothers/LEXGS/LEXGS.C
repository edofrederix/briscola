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
    symmetric_(dict.lookupOrDefault<Switch>("symmetric", true))
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
    xi = singular ? gAverage(x) : List<Type>(x.size(), Zero);

    const List<bool> diagonal(sys.diagonal());

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
                meshDirection<Type,MeshType>& xd = x[d];

                const meshDirection<SType,MeshType>& Ad = A[d];
                const meshDirection<Type,MeshType>& bd = b[d];

                forAllCells(xd, i, j, k)
                {
                    Switch forcing = false;

                    forAll(sys.x().immersedBoundaryConditions(), ib)
                    {
                        if
                        (
                            sys.x().immersedBoundaryConditions()[ib]
                                .forcingPoints()(l,d,i,j,k)
                        )
                        {
                            forcing = true;
                        }
                    }

                    if (!forcing)
                    {
                        xd(i,j,k) =
                            (
                            bd(i,j,k)
                            - lowerRowProduct(Ad,xd,i,j,k)
                            - upperRowProduct(Ad,xd,i,j,k)
                            - xi[d]
                            )
                        / Ad(i,j,k).center();
                    }
                }

                    if (!symmetric_)
                        continue;

                forAllCells(xd, i, j, k)
                {
                    Switch forcing = false;

                    forAll(sys.x().immersedBoundaryConditions(), ib)
                    {
                        if
                        (
                            sys.x().immersedBoundaryConditions()[ib]
                                .forcingPoints()(l,d,i,j,k)
                        )
                        {
                            forcing = true;
                        }
                    }

                    if (!forcing)
                    {
                        xd(i,j,k) =
                            (
                            bd(i,j,k)
                            - lowerRowProduct(Ad,xd,i,j,k)
                            - upperRowProduct(Ad,xd,i,j,k)
                            - xi[d]
                            )
                        / Ad(i,j,k).center();
                    }
                }
            }
        }

        x.correctNonEliminatedBoundaryConditions();
        if (l == 0)
            x.correctImmersedBoundaryConditions();
    }

    x.correctEliminatedBoundaryConditions();
    x.correctUnsetBoundaryConditions();

    if (!singular)
        xi.clear();
}

}

}

}
