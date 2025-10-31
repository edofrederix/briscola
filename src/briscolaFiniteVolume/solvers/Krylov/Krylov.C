#include "Krylov.H"
#include "diagonal.H"
#include "diagonalSmoother.H"
#include "PstreamGlobals.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
void Krylov<SType,Type,MeshType>::prepare
(
    linearSystem<SType,Type,MeshType>& sys
)
{
    const labelVector* offsets = SType::offsets;

    const MPI_Comm& comm =
        Pstream::parRun()
      ? PstreamGlobals::MPI_COMM_FOAM
      : PETSC_COMM_SELF;

    const meshLevel<SType,MeshType>& A = sys.A()[0];

    const meshLevel<label,MeshType>& numbers =
        sys.fvMsh().template metrics<MeshType>().globalCellNumbers()[0];

    const List<bool> diagonal(sys.diagonal());

    forAll(solvers_, i)
        if (solvers_.set(i))
            KSPDestroy(&solvers_[i]);

    forAll(matrices_, i)
        if (matrices_.set(i))
            MatDestroy(&matrices_[i]);

    forAll(A, d)
    if (!diagonal[d])
    {
        const label m = A[d].size();
        const label M = returnReduce(m, sumOp<label>());
        const label nz = Foam::min(label(SType::nCsComponents),M);

        // Create matrix

        matrices_.set(d, new Mat);
        Mat& mat = matrices_[d];

        MatCreate(comm, &mat);
        MatSetSizes(mat, m, m, M, M);

        if (Pstream::parRun())
        {
            MatSetType(mat, MATMPIAIJ);
            MatMPIAIJSetPreallocation(mat, nz, NULL, nz, NULL);
        }
        else
        {
            MatSetType(mat, MATSEQAIJ);
            MatSeqAIJSetPreallocation(mat, nz, NULL);
        }

        MatSetUp(mat);

        forAllCells(A[d], i, j, k)
        {
            const labelVector ijk(i,j,k);

            SType coeffs = A[d](ijk);

            labelList colNums(SType::nCsComponents);

            forAll(colNums, s)
                colNums[s] = numbers(d, ijk+offsets[s]);

            // Remove boundary coefficients

            forAll(colNums, s)
                if (colNums[s] < 0 || colNums[s] >= M)
                    coeffs[s] = 0;

            // Move redundant indices

            forAll(colNums, s)
            {
                const int t = findIndex(colNums, colNums[s]);

                if (t < s)
                {
                    coeffs[t] += coeffs[s];
                    coeffs[s] = 0;
                }
            }

            // Add to matrix

            forAll(colNums, s)
                if (coeffs[s] != 0)
                    MatSetValue
                    (
                        mat,
                        colNums[0],
                        colNums[s],
                        coeffs[s],
                        INSERT_VALUES
                    );
        }

        MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);

        // Setup solver

        solvers_.set(d, new KSP);
        KSP& ksp = solvers_[d];

        KSPCreate(comm, &ksp);
        KSPSetOperators(ksp, mat, mat);
        KSPSetType(ksp, solverType_.c_str());
        KSPSetNormType(ksp, KSP_NORM_UNPRECONDITIONED);
    }

    this->prepared_ = true;
}

template<class SType, class Type, class MeshType>
void Krylov<SType,Type,MeshType>::solve
(
    linearSystem<SType,Type,MeshType>& sys,
    const scalar relTol,
    const scalar absTol,
    const label maxIter
)
{
    const labelVector* offsets = SType::offsets;
    const label nDir = MeshType::numberOfDirections;

    meshLevel<Type, MeshType>& x = sys.x()[0];
    const meshLevel<Type, MeshType>& b = sys.b()[0];
    const meshLevel<SType, MeshType>& A = sys.A()[0];

    const meshLevel<label,MeshType>& numbers =
        sys.fvMsh().template metrics<MeshType>().globalCellNumbers()[0];

    const MPI_Comm& comm =
        Pstream::parRun()
      ? PstreamGlobals::MPI_COMM_FOAM
      : PETSC_COMM_SELF;

    const label n = list(pTraits<Type>::one).size();

    // Correct the boundary conditions

    x.correctBoundaryConditions();

    // Buffer

    meshLevel<Type, MeshType> buffer(sys.fvMsh(), 0);

    // Initial residual

    sys.residual(buffer);

    // Residual normalization factors

    const List<Type> normFactors(this->normFactors(sys,buffer));

    const List<Type> initialResiduals =
        cmptDivide
        (
            Foam::cmptSqrt(gSum(cmptSqr(buffer))),
            normFactors
        );

    const List<Type> bNormFactors
    (
        Foam::cmptSqrt(gSum(cmptSqr(b)))
    );

    for (int d = 0; d < nDir; d++)
    {
        const label m = A[d].size();
        const label M = returnReduce(m, sumOp<label>());

        Vec sol, rhs;

        VecCreateMPI(comm, m, M, &sol);
        VecCreateMPI(comm, m, M, &rhs);

        VecSetUp(sol);
        VecSetUp(rhs);

        // Copy to buffer

        buffer[d] = b[d];

        // Add boundary sources

        forAllCells(A, i, j, k)
        {
            const labelVector ijk(i,j,k);

            const SType coeffs = A[d](ijk);

            labelList colNums(SType::nCsComponents);

            forAll(colNums, s)
                colNums[s] = numbers[d](ijk + offsets[s]);

            forAll(colNums, s)
                if (colNums[s] < 0 || colNums[s] >= M)
                    buffer(d,i,j,k) -= x(d,i,j,k)*coeffs[s];
        }

        // We need to solve each vector direction separately

        KSP& ksp = solvers_[d];

        for (int dir = 0; dir < n; dir++)
        {
            // Copy buffer to PETSc right-hand side type

            label row = numbers(d,0,0,0);
            forAllCells(b, i, j, k)
            {
                VecSetValue
                (
                    rhs,
                    row,
                    scalar_cast(&buffer(d,i,j,k))[dir],
                    INSERT_VALUES
                );

                VecSetValue
                (
                    sol,
                    row,
                    scalar_cast(&x(d,i,j,k))[dir],
                    INSERT_VALUES
                );

                row++;
            }

            VecAssemblyBegin(rhs);
            VecAssemblyEnd(rhs);

            VecAssemblyBegin(sol);
            VecAssemblyEnd(sol);

            // Set tolerances

            const scalar nf = scalar_cast(&normFactors[d])[dir];
            const scalar bnf = scalar_cast(&bNormFactors[d])[dir];
            const scalar initialResidual =
                scalar_cast(&initialResiduals[d])[dir];

            KSPSetTolerances
            (
                ksp,
                relTol*nf*initialResidual/Foam::max(bnf, 1e-12),
                absTol*nf,
                1e3,
                maxIter
            );

            // Solve

            KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
            KSPSolve(ksp, rhs, sol);

            label nIter;
            KSPGetIterationNumber(ksp, &nIter);

            KSPConvergedReason reason;
            KSPGetConvergedReason(ksp, &reason);

            if (reason < 0)
                WarningInFunction
                    << "PETSc solver"
                    << " did not converge with "
                    << "reason " << reason << endl;

            scalar finalResidual;
            KSPGetResidualNorm(ksp, &finalResidual);

            finalResidual /= nf;

            // Copy back to x

            PetscScalar *arr;
            VecGetArray(sol, &arr);

            row = 0;
            forAllCells(x[d], i, j, k)
                scalar_cast(&x(d,i,j,k))[dir] = arr[row++];

            VecRestoreArray(sol, &arr);

            Info<< this->typeName
                << ": Solving for " << MeshType::typeName
                << " " << sys.x().name()
                << (nDir > 1 ? word(" direction " + Foam::name(d)) : word(""))
                << (n > 1 ? word(" component " + Foam::name(dir)) : word(""))
                << ", initial residual = " << initialResidual
                << ", final residual = " << finalResidual
                << ", nIter = " << nIter << endl;
        }

        VecDestroy(&sol);
        VecDestroy(&rhs);
    }

    x.correctBoundaryConditions();
}

template<class SType, class Type, class MeshType>
Krylov<SType,Type,MeshType>::Krylov
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    solver<SType,Type,MeshType>(dict,fvMsh),
    solverType_(dict.lookupOrDefault<word>("solverType", KSPBCGS)),
    matrices_(MeshType::numberOfDirections),
    solvers_(MeshType::numberOfDirections),
    prepared_(false)
{
    // Initialize PETSc

    if (Pstream::parRun())
        PETSC_COMM_WORLD = PstreamGlobals::MPI_COMM_FOAM;

    PetscBool initialized;
    PetscInitialized(&initialized);

    if (!initialized)
        PetscInitialize(NULL, NULL, NULL, NULL);
}

template<class SType, class Type, class MeshType>
Krylov<SType,Type,MeshType>::~Krylov()
{}

template<class SType, class Type, class MeshType>
void Krylov<SType,Type,MeshType>::solve
(
    linearSystem<SType,Type,MeshType>& sys,
    const bool constMatrix
)
{
    const label nDir = MeshType::numberOfDirections;

    if (SType::nCsComponents > 1)
        sys.eliminateGhosts();

    sys.setForcingMask();

    if
    (
        SType::nCsComponents == 1
     || sum(sys.diagonal()) == nDir
    )
    {
        sys.x().makeShallow();
        sys.b().makeShallow();

        diagonalSmoother<SType,Type,MeshType>::Smooth
        (
            sys,
            0,
            1,
            labelList(MeshType::numberOfDirections, 0)
        );

        this->printSolverStats
        (
            diagonal<SType,Type,MeshType>::typeName,
            sys.x().name(),
            List<Type>(nDir, pTraits<Type>::one),
            List<Type>(nDir, Zero),
            0
        );
    }
    else
    {
        if (!constMatrix || !this->prepared_)
            this->prepare(sys);

        this->solve
        (
            sys,
            this->relTol_,
            this->tolerance_,
            this->maxIter_
        );
    }
}

}

}

}
