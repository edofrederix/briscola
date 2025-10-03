#include "PETScExternalSolver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

PETScExternalSolver::PETScExternalSolver(const dictionary& dict)
:
    PETScSolverBase(dict)
{}

PETScExternalSolver::PETScExternalSolver(const PETScExternalSolver& s)
:
    PETScSolverBase(s)
{}

PETScExternalSolver::~PETScExternalSolver()
{}

void PETScExternalSolver::prepare(const Mat& mat, const MPI_Comm& comm)
{
    if (factor_.valid())
        MatDestroy(&factor_());

    factor_.set(new Mat);
}

void PETScExternalSolver::solve
(
    Vec& x,
    const Vec& rhs,
    const MPI_Comm& comm
)
{
    if (&x == &rhs)
    {
        Vec tmp;
        VecDuplicate(rhs, &tmp);

        MatSolve(factor_(), rhs, tmp);

        VecCopy(tmp, x);
        VecDestroy(&tmp);
    }
    else
    {
        MatSolve(factor_(), rhs, x);
    }
}

}

}

}
