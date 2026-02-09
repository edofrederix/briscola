#include "APLU.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
APLU<SType,Type,MeshType>::APLU
(
    const dictionary& dict,
    const fvMesh& fvMsh,
    const label l
)
:
    solver<SType,Type,MeshType>::externalSolver(dict,fvMsh,l),
    APtr_(),
    DPtrs_(MeshType::numberOfDirections),
    pivotPtrs_(MeshType::numberOfDirections)
{}

template<class SType, class Type, class MeshType>
void APLU<SType,Type,MeshType>::prepare(linearSystem<SType,Type,MeshType>& sys)
{
    if (!APtr_.valid())
    {
        APtr_.reset
        (
            new linearSystemAggregation<SType,Type,MeshType>
            (
                sys,
                this->l_,
                1
            )
        );
    }

    const linearSystemAggregation<SType,Type,MeshType>& lsa = APtr_();

    if (!sys.diagonal())
    {
        forAll(pivotPtrs_, d)
        {
            // Get full matrix in compressed row format

            scalarList values;
            labelList inners;
            labelList outers;

            lsa.compressedRowFormat
            (
                values,
                inners,
                outers,
                sys,
                d,
                true
            );

            if (lsa.master())
            {
                const label n = outers.size()-1;

                // Create dense matrix

                DPtrs_.set(d, new scalarSquareMatrix(n));
                scalarSquareMatrix& D = DPtrs_[d];
                D = Zero;

                // Fill dense matrix

                for (int i = 0; i < outers.size()-1; i++)
                    for (int j = outers[i]; j < outers[i+1]; j++)
                        D[i][inners[j]] = values[j];

                pivotPtrs_.set(d, new labelList(n));
                labelList& pivots = pivotPtrs_[d];

                LUDecompose(D, pivots);
            }
        }
    }

    solver<SType,Type,MeshType>::externalSolver::prepare(sys);
}

template<class SType, class Type, class MeshType>
void APLU<SType,Type,MeshType>::solve(linearSystem<SType,Type,MeshType>& sys)
{
    if (!this->prepared_)
        this->prepare(sys);

    const linearSystemAggregation<SType,Type,MeshType>& lsa = APtr_();

    meshLevel<Type,MeshType>& x = sys.x()[this->l_];
    const meshLevel<Type,MeshType>& b = sys.b()[this->l_];
    const meshLevel<SType,MeshType>& A = sys.A()[this->l_];

    for (int d = 0; d < MeshType::numberOfDirections; d++)
    {
        if (sys.diagonal())
        {
            forAllCells(x[d], i, j, k)
                x[d](i,j,k) = b[d](i,j,k)/A[d](i,j,k).center();
        }
        else
        {
            List<Type> buffer;
            lsa.rhsSource(buffer, sys, d);

            if (lsa.master())
            {
                scalarSquareMatrix& D = DPtrs_[d];
                labelList& pivots = pivotPtrs_[d];

                LUBacksubstitute(D, pivots, buffer);
            }

            lsa.distribute(x[d], buffer);
        }
    }

    x.correctBoundaryConditions();
}

}

}

}
