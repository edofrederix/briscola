#include "PETScDirectSolver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

PETScDirectSolver::PETScDirectSolver(const dictionary& dict)
:
    PETScSolverBase(dict)
{}

PETScDirectSolver::PETScDirectSolver(const PETScDirectSolver& s)
:
    PETScSolverBase(s)
{}

PETScDirectSolver::~PETScDirectSolver()
{}

void PETScDirectSolver::prepare(const Mat& mat, const MPI_Comm& comm)
{
    if (factor_.valid())
        MatDestroy(factor_.ptr());

    factor_.set(new Mat);
}

void PETScDirectSolver::solve
(
    Vec& x,
    const Vec& rhs,
    const MPI_Comm& comm
)
{
    if (&x == &rhs)
    {
        Vec x;
        VecDuplicate(rhs, &x);

        MatSolve(factor_(), rhs, x);

        VecCopy(x, rhs);
    }
    else
    {
        MatSolve(factor_(), rhs, x);
    }
}

}

}

}
