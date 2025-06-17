#include "PETScKSPSolver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

PETScKSPSolver::PETScKSPSolver(const dictionary& dict)
:
    PETScSolverBase(dict),
    tolerance_(dict.lookupOrDefault<scalar>("tolerance", 0)),
    relTol_(dict.lookupOrDefault<scalar>("relTol", 1e-5)),
    maxIter_(dict.lookupOrDefault<label>("maxIter", 1000)),
    printStats_(dict.lookupOrDefault<Switch>("printStats", false))
{}

PETScKSPSolver::PETScKSPSolver(const PETScKSPSolver& s)
:
    PETScSolverBase(s),
    tolerance_(s.tolerance_),
    relTol_(s.relTol_),
    maxIter_(s.maxIter_),
    printStats_(s.printStats_)
{}

PETScKSPSolver::~PETScKSPSolver()
{}

void PETScKSPSolver::prepare(const Mat& mat, const MPI_Comm& comm)
{
    solverPtr_.reset(new KSP);
    KSP& ksp = solverPtr_();

    KSPCreate(comm, &ksp);
    KSPSetOperators(ksp, mat, mat);
    KSPSetTolerances(ksp, relTol_, tolerance_, 1e3, maxIter_);
}

void PETScKSPSolver::solve
(
    Vec& x,
    const Vec& rhs,
    const MPI_Comm& comm
)
{
    KSP& ksp = solverPtr_();

    KSPSolve(ksp, x, rhs);

    KSPConvergedReason reason;
    KSPGetConvergedReason(ksp, &reason);

    if (reason < 0)
    {
        WarningInFunction
            << "PETSc solver " << this->type() << " did not converge with "
            << "reason " << reason << endl;
    }
    if (printStats_)
    {
        PetscInt nIter;
        KSPGetIterationNumber(ksp, &nIter);

        scalar finalResidual;
        KSPGetResidualNorm(ksp, &finalResidual);

        // Initial guess is zero so the initial residual is unity

        Info<< "PETSc " << this->type() << " solver"
            << ": Initial residual = 1"
            << ", final residual = " << finalResidual
            << ", convergence reason = " << reason
            << ", nIter = " << nIter << endl;
    }
}

}

}

}
