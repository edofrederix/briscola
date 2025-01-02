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
    nAggregationParts_(dict.lookupOrDefault<label>("nAggregationParts", 1)),
    minIter_(dict.lookupOrDefault<label>("minIter", 1)),
    maxIter_
    (
        nAggregationParts_ == 1
      ? 1
      : dict.lookupOrDefault<label>("maxIter", nAggregationParts_)
    ),
    relTol_(dict.lookupOrDefault<scalar>("relTol", 1e-3)),
    tolerance_(dict.lookupOrDefault<scalar>("tolerance", 0)),
    printStats_(dict.lookupOrDefault<Switch>("printStats", false))
{
    if (nAggregationParts_ > Pstream::nProcs())
        FatalErrorInFunction
            << "Requested number of aggregation parts is larger "
            << "than the number of processors" << endl
            << abort(FatalError);

    for (int d = 0; d < MeshType::numberOfDirections; d++)
        solverPtrs_.set
        (
            d,
            EigenSolver::New
            (
                dict.lookupOrDefault<word>("EigenSolver", "default")
            ).ptr()
        );
}

template<class SType, class Type, class MeshType>
void Eigen<SType,Type,MeshType>::prepare
(
    linearSystem<SType,Type,MeshType>& sys
)
{
    if (APtr_.valid())
    {
        // The linear system is already initialized so update the matrix only

        APtr_->update(sys);
    }
    else
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

    meshLevel<Type,MeshType> r(sys.fvMsh(), this->l_);
    sys.residual(r);

    normFactors_ =
        solver<SType,Type,MeshType>::normFactors(sys,r);

    initialResiduals_ =
        cmptDivide(gSum(cmptMag(r)), normFactors_);

    const List<bool> diagonal(sys.diagonal());

    if (APtr_->master())
        forAll(solverPtrs_, d)
            if (!diagonal[d])
                solverPtrs_[d].prepare
                (
                    APtr_->operator[](d),
                    Foam::max
                    (
                        0.5*Foam::cmptMax(initialResiduals_[d])*relTol_,
                        0.5*tolerance_
                    )
                );

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

    const List<bool> singular(sys.singular());
    const List<bool> diagonal(sys.diagonal());

    meshLevel<Type,MeshType>& x = sys.x()[this->l_];
    const meshLevel<Type,MeshType>& b = sys.b()[this->l_];
    const meshLevel<SType,MeshType>& A = sys.A()[this->l_];

    List<Type> currentResiduals(initialResiduals_);

    labelList converged =
        solver<SType,Type,MeshType>::checkConvergence
        (
            currentResiduals,
            initialResiduals_,
            relTol_,
            tolerance_
        );

    label iter = 0;

    // Set right-hand sides

    List<List<Type>> rhs(MeshType::numberOfDirections);

    for (int d = 0; d < MeshType::numberOfDirections; d++)
        if (!diagonal[d])
            E.rhsSource(rhs[d], A[d], x[d], b[d]);

    List<List<Type>> solutions(MeshType::numberOfDirections);

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
                    // Copy rhs to Eigen right-hand side type

                    EigenSolver::EigenRhsType B(rhs[d].size(),n);

                    forAll(rhs[d], i)
                        for (int j = 0; j < n; j++)
                            B(i,j) = scalar_cast(&rhs[d][i])[j];

                    // Solve

                    EigenSolver::EigenRhsType X(rhs[d].size(),n);

                    if (solverPtrs_[d].needGuess())
                        forAll(solutions[d], i)
                            for (int j = 0; j < n; j++)
                                X(i,j) =
                                    solutions[d].size()
                                  ? scalar_cast(&solutions[d][i])[j]
                                  : 0;

                    solverPtrs_[d].solve(X,B);

                    // Copy back to buffer

                    solutions[d].resize(rhs[d].size());

                    forAll(solutions[d], i)
                        for (int j = 0; j < n; j++)
                            scalar_cast(&solutions[d][i])[j] = X(i,j);

                    // Send/copy solutions

                    label offset = 0;
                    forAll(E.colNums()[d], proc)
                    {
                        const label procNum = E.myPartMasterNum() + proc;

                        if (Pstream::myProcNo() == procNum)
                        {
                            label c = 0;
                            forAllCells(x[d], i, j, k)
                                x[d](i,j,k) = solutions[d][offset + c++];
                        }
                        else
                        {
                            UOPstream::write
                            (
                                Pstream::commsTypes::blocking,
                                procNum,
                                reinterpret_cast<char*>(&solutions[d][offset]),
                                E.colNums()[d][proc].size()*sizeof(Type),
                                0,
                                UPstream::worldComm
                            );
                        }

                        offset += E.colNums()[d][proc].size();
                    }
                }
                else
                {
                    // Receive and copy solution

                    solutions[d].resize(x[d].size());

                    UIPstream::read
                    (
                        Pstream::commsTypes::blocking,
                        E.myPartMasterNum(),
                        reinterpret_cast<char*>(solutions[d].begin()),
                        solutions[d].byteSize(),
                        0,
                        UPstream::worldComm
                    );

                    label c = 0;
                    forAllCells(x[d], i, j, k)
                        x[d](i,j,k) = solutions[d][c++];
                }
            }
        }

        // Assuming that the solution for each aggregation matrix is exact, we
        // can recompute the residual from the change in the right-hand side
        // update after the aggregation matrix boundary condition correction

        E.correctBoundaryConditions(x);

        forAll(x, d)
        {
            List<Type> r(rhs[d]);

            if (!diagonal[d] && !converged[d])
            {
                E.rhsSource(rhs[d], A[d], x[d], b[d]);

                forAll(r, i)
                    r[i] -= rhs[d][i];

                currentResiduals[d] =
                    Foam::cmptDivide(gSum(cmptMag(r)), normFactors_[d]);
            }
        }

        converged =
            solver<SType,Type,MeshType>::checkConvergence
            (
                currentResiduals,
                initialResiduals_,
                relTol_,
                tolerance_
            );

        iter++;
    }

    // Set average to zero in singular directions

    if (sum(singular) > 0)
        for (int d = 0; d < MeshType::numberOfDirections; d++)
            if (singular[d])
                x[d] -= gAverage(x[d]);

    x.correctBoundaryConditions();

    // Print solver statistics

    if (printStats_)
        solver<SType,Type,MeshType>::printSolverStats
        (
            this->type() + " " + solverPtrs_[0].type(),
            sys.x().name(),
            initialResiduals_,
            currentResiduals,
            iter
        );
}

}

}

}
