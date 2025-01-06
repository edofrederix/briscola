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
    solverPtrs_(MeshType::numberOfDirections),
    nAggregationParts_
    (
        dict.lookupOrDefault<label>
        (
            "nAggregationParts",
            1
        )
    ),
    minIter_(dict.lookupOrDefault<label>("minIter", 1)),
    maxIter_
    (
        nAggregationParts_ == 1
      ? 1
      : dict.lookupOrDefault<label>("maxIter", 100)
    ),
    relTol_(dict.lookupOrDefault<scalar>("relTol", 1e-4)),
    tolerance_(dict.lookupOrDefault<scalar>("tolerance", 0)),
    printStats_(dict.lookupOrDefault<Switch>("printStats", false))
{
    if (nAggregationParts_ > Pstream::nProcs())
        FatalErrorInFunction
            << "Requested number of aggregation parts is larger "
            << "than the number of processors" << endl
            << abort(FatalError);

    forAll(solverPtrs_, d)
        solverPtrs_.set
        (
            d,
            EigenSolverBase::New
            (
                dict.lookupOrDefault<word>("EigenSolver", "SparseLU"),
                dict.lookupOrDefault<dictionary>
                (
                    "EigenSolverCoeffs",
                    dictionary::null
                )
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
                nAggregationParts_
            )
        );
    }

    const List<bool> diagonal(sys.diagonal());

    forAll(solverPtrs_, d)
    if (!diagonal[d])
    {
        APtr_->update(solverPtrs_[d].order(), d);

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

    // Collect right-hand sides, residuals and initial solutions

    List<List<Type>> rhs(MeshType::numberOfDirections);
    List<List<Type>> sol(MeshType::numberOfDirections);

    for (int d = 0; d < MeshType::numberOfDirections; d++)
    {
        if (!diagonal[d])
        {
            E.rhsSource(rhs[d], sys, d);
            E.collect(sol[d], x[d]);
        }
    }

    // Residual normalization factors, assuming the initial solution is zero.
    // All processors must take part in this.

    List<Type> normFactors
    (
        MeshType::numberOfDirections,
        1e-20*pTraits<Type>::one
    );

    for (int d = 0; d < MeshType::numberOfDirections; d++)
    {
        if (!diagonal[d] && (maxIter_ > 1 || printStats_))
        {
            Type f = Zero;

            forAllCells(b[d], i, j, k)
                f += Foam::cmptSqr(b[d](i,j,k));

            f = Foam::cmptSqrt(returnReduce(f, sumOp<Type>()));
            normFactors[d] = Foam::max(f, normFactors[d]);
        }
    }

    const List<Type> initialResiduals
    (
        MeshType::numberOfDirections,
        pTraits<Type>::one
    );

    List<Type> currentResiduals(initialResiduals);
    labelList converged(MeshType::numberOfDirections, 0);

    label iter = 0;

    List<List<Type>> prev(MeshType::numberOfDirections);

    while
    (
        (iter < maxIter_ && Foam::min(converged) == 0)
     || iter < minIter_
    )
    {
        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            if (converged[d])
            {
                // Nothing to do
            }
            if (diagonal[d])
            {
                forAllCells(x[d], i, j, k)
                    x[d](i,j,k) = b[d](i,j,k)/A[d](i,j,k).center();
            }
            else
            {
                if (E.master())
                {
                    // Cache the solution for later use

                    if (iter < maxIter_-1 || printStats_)
                        prev[d] = sol[d];

                    // Copy rhs to Eigen right-hand side type

                    EigenSolverBase::RhsType B(rhs[d].size(),n);

                    forAll(rhs[d], i)
                        for (int j = 0; j < n; j++)
                            B(i,j) = scalar_cast(&rhs[d][i])[j];

                    // Solve

                    EigenSolverBase::RhsType X(rhs[d].size(),n);

                    if (solverPtrs_[d].needGuess())
                        forAll(sol[d], i)
                            for (int j = 0; j < n; j++)
                                X(i,j) = scalar_cast(&sol[d][i])[j];

                    solverPtrs_[d].solve(X,B);

                    // Copy back to buffer

                    sol[d].resize(rhs[d].size());

                    forAll(sol[d], i)
                        for (int j = 0; j < n; j++)
                            scalar_cast(&sol[d][i])[j] = X(i,j);
                }

                E.distribute(x[d], sol[d]);
            }
        }

        E.correctBoundaryConditions(x);

        forAll(x, d)
        {
            if
            (
                !diagonal[d]
             && !converged[d]
             && (iter < maxIter_-1 || printStats_)
            )
            {
                // Update right-hand side with updated boundary sources

                E.rhsSource(rhs[d], sys, d);

                // The residual is determined from the change in the solution

                Type f = Zero;
                forAll(sol[d], i)
                    f += Foam::cmptSqr(sol[d][i] - prev[d][i]);

                currentResiduals[d] =
                    Foam::cmptDivide
                    (
                        Foam::cmptSqrt(returnReduce(f, sumOp<Type>())),
                        normFactors[d]
                    );
            }
        }

        converged =
            solver<SType,Type,MeshType>::checkConvergence
            (
                currentResiduals,
                initialResiduals,
                relTol_,
                tolerance_
            );

        iter++;
    }

    x.correctBoundaryConditions();

    // Print solver statistics

    if (printStats_)
        solver<SType,Type,MeshType>::printSolverStats
        (
            this->type() + " " + solverPtrs_[0].type(),
            sys.x().name(),
            initialResiduals,
            currentResiduals,
            iter
        );
}

}

}

}
