#include "Eigen.H"

#include "SquareMatrix.H"
#include "LUscalarMatrix.H"
#include "meshDirectionStencilFunctions.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
Eigen<SType,Type,MeshType>::Eigen
(
    const dictionary& dict,
    const fvMesh& fvMsh,
    const label l
)
:
    solver<SType,Type,MeshType>::directSolver(dict,fvMsh,l),
    APtr_(),
    solverPtrs_(MeshType::numberOfDirections)
{
    forAll(solverPtrs_, d)
        solverPtrs_.set
        (
            d,
            EigenSolverBase::New
            (
                dict.lookupOrDefault<word>("EigenSolver", "SparseLU"),
                dict
            ).ptr()
        );
}

template<class SType, class Type, class MeshType>
void Eigen<SType,Type,MeshType>::prepare
(
    linearSystem<SType,Type,MeshType>& sys
)
{
    if (!APtr_.valid())
    {
        APtr_.reset
        (
            new EigenLinearSystem<SType,Type,MeshType>
            (
                sys,
                this->l_,
                1
            )
        );
    }

    const List<bool> diagonal(sys.diagonal());

    forAll(solverPtrs_, d)
    {
        if (!diagonal[d])
        {
            APtr_->prepare(solverPtrs_[d].order(), d);

            if (APtr_->master())
            {
                if (solverPtrs_[d].order() == ::Eigen::RowMajor)
                {
                    solverPtrs_[d].prepare(APtr_->rowMajorMatrix(d));
                }
                else
                {
                    solverPtrs_[d].prepare(APtr_->colMajorMatrix(d));
                }
            }
        }
    }

    solver<SType,Type,MeshType>::directSolver::prepare(sys);
}

template<class SType, class Type, class MeshType>
void Eigen<SType,Type,MeshType>::solve
(
    linearSystem<SType,Type,MeshType>& sys
)
{
    if (!this->prepared_)
        this->prepare(sys);

    const EigenLinearSystem<SType,Type,MeshType>& E = APtr_();

    const label n = list(pTraits<Type>::one).size();

    const List<bool> diagonal(sys.diagonal());

    meshLevel<Type,MeshType>& x = sys.x()[this->l_];
    const meshLevel<Type,MeshType>& b = sys.b()[this->l_];
    const meshLevel<SType,MeshType>& A = sys.A()[this->l_];

    for (int d = 0; d < MeshType::numberOfDirections; d++)
    {
        if (diagonal[d])
        {
            forAllCells(x[d], i, j, k)
                x[d](i,j,k) = b[d](i,j,k)/A[d](i,j,k).center();
        }
        else
        {
            List<Type> buffer;
            E.rhsSource(buffer, sys, d);

            if (E.master())
            {
                // Copy buffer to Eigen right-hand side type

                EigenSolverBase::RhsType B(buffer.size(),n);

                forAll(buffer, i)
                    for (int j = 0; j < n; j++)
                        B(i,j) = scalar_cast(&buffer[i])[j];

                // Solve

                EigenSolverBase::RhsType X;
                solverPtrs_[d].solve(X,B);

                // Copy back to buffer

                forAll(buffer, i)
                    for (int j = 0; j < n; j++)
                        scalar_cast(&buffer[i])[j] = X(i,j);
            }

            E.distribute(x[d], buffer);
        }
    }

    x.correctBoundaryConditions();
}

}

}

}
