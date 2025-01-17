#include "PETSc.H"
#include "PstreamGlobalsLsa.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
PETSc<SType,Type,MeshType>::PETSc
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
    )
{
    // Initialize PETSc

    if (Pstream::parRun())
        PETSC_COMM_WORLD = PstreamGlobals::MPI_COMM_FOAM;

    PetscBool initialized;
    PetscInitialized(&initialized);

    if (!initialized)
        PetscInitialize(NULL, NULL, NULL, NULL);

    // Setup solvers

    forAll(solverPtrs_, d)
        solverPtrs_.set
        (
            d,
            PETScSolverBase::New
            (
                dict.lookupOrDefault<word>("PETScSolver", "PCLU"),
                dict
            ).ptr()
        );

    // Force aggregation parts to one if the solver is not parallel

    if (!solverPtrs_[0].parallel())
        nAggregationParts_ = 1;
}

template<class SType, class Type, class MeshType>
PETSc<SType,Type,MeshType>::~PETSc()
{}

template<class SType, class Type, class MeshType>
void PETSc<SType,Type,MeshType>::prepare
(
    linearSystem<SType,Type,MeshType>& sys
)
{
    if (!APtr_.valid())
    {
        APtr_.reset
        (
            new PETScLinearSystem<SType,Type,MeshType>
            (
                sys,
                this->l_,
                nAggregationParts_
            )
        );
    }

    const List<bool> diagonal(sys.diagonal());

    const MPI_Comm& comm =
        Pstream::parRun()
      ? PstreamGlobals::lsaGetComm(APtr_->masterCommNum())
      : PETSC_COMM_SELF;

    forAll(solverPtrs_, d)
    {
        if (!diagonal[d])
        {
            APtr_->prepare(d);

            if (APtr_->master())
                solverPtrs_[d].prepare
                (
                    APtr_->matrix(d),
                    comm
                );
        }
    }

    solver<SType,Type,MeshType>::directSolver::prepare(sys);
}

template<class SType, class Type, class MeshType>
void PETSc<SType,Type,MeshType>::solve
(
    linearSystem<SType,Type,MeshType>& sys
)
{
    if (!this->prepared_)
        this->prepare(sys);

    const PETScLinearSystem<SType,Type,MeshType>& P = APtr_();

    const MPI_Comm& comm =
        Pstream::parRun()
      ? PstreamGlobals::lsaGetComm(P.masterCommNum())
      : PETSC_COMM_SELF;

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
            P.rhsSource(buffer, sys, d);

            if (P.master())
            {
                Vec rhs;

                VecCreateMPI
                (
                    comm,
                    P.partSize(d),
                    P.globalSize(d),
                    &rhs
                );

                VecSetUp(rhs);

                // We need to solve each vector direction separately

                for (int dir = 0; dir < n; dir++)
                {
                    // Copy buffer to PETSc right-hand side type

                    label row = P.partStart(d);
                    forAll(buffer, i)
                        VecSetValue
                        (
                            rhs,
                            row++,
                            scalar_cast(&buffer[i])[dir],
                            INSERT_VALUES
                        );

                    VecAssemblyBegin(rhs);
                    VecAssemblyEnd(rhs);

                    // Solve

                    solverPtrs_[d].solve
                    (
                        rhs,
                        rhs,
                        comm
                    );

                    // Copy back to buffer

                    PetscScalar *arr;
                    VecGetArray(rhs, &arr);

                    row = 0;
                    forAll(buffer, i)
                        scalar_cast(&buffer[i])[dir] = arr[row++];

                    VecRestoreArray(rhs, &arr);
                }

                VecDestroy(&rhs);
            }

            P.distribute(x[d], buffer);
        }
    }

    x.correctBoundaryConditions();
}

}

}

}
